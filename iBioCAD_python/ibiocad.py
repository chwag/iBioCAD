######!/usr/bin/python

print("Content-type: text/html\n\n")

'''
Add modules uploaded to the cPanel file manager to the system path.
Will need to be modified depending on the file manager/server.
May not be needed depending on the nature of the file manager/server
(isn't needed if developer has access to usr/bin/python)
'''
import sys
#sys.path.append('/home/ibiocadsite/public_html/cgi-bin/iBioCAD/lib/python2.6/site-packages')
try:
    reload(sys)
    sys.setdefaultencoding('UTF8')
except:
    pass

'''
Create a virtual envrionment for the project if desired.
Allows the project to run withoutthe dependencies of the system.
venv is standard in python3, can be backported to python2
'''
###import venv
###venv.create('/home/ibiocadsite/public_html/cgi-bin/iBioCAD')


import webapp2      #web framework
import os
import jinja2       #optional template rendering
from Bio.SeqUtils import MeltingTemp as mt
from Bio import SeqIO
is_local = False     #set file path to the local folder if True or server path if False
import itertools
import io

#import stylesheets and javascript for web pages
if is_local:
    #path on local
    with open('templates/main_page.css','r') as my_main_page_css_raw:
        my_main_page_css = my_main_page_css_raw.read()
    with open('templates/semantic.min.css','r') as semantic_css_raw:
        semantic_css = semantic_css_raw.read()
    with open('templates/semantic.min.js','r') as semantic_js_raw:
        semantic_js = semantic_js_raw.read()
    css = my_main_page_css + semantic_css
    js = semantic_js
else:
    #path on the server
    with open('/var/www/ibiocad/iBioCAD/templates/main_page.css','r') as my_main_page_css_raw:
        my_main_page_css = my_main_page_css_raw.read()
    with open('/var/www/ibiocad/iBioCAD/templates/semantic.min.css','r') as semantic_css_raw:
        semantic_css = semantic_css_raw.read()
    with open('/var/www/ibiocad/iBioCAD/templates/semantic.min.js','r') as semantic_js_raw:
        semantic_js = semantic_js_raw.read()
    css = my_main_page_css + semantic_css
    js = semantic_js

tkinter_use = False     #set to true to use a GUI window within a web handler
server_debug = False     #set to false to disable the ability to restart a server process
'''
Initialize jinja2 rendering
   -html documents that will be rendered must be put in a file folder named "templates"
   -jinja2 templates add functionality with special syntax.
'''
template_dir = os.path.join(os.path.dirname(__file__), "templates")
jinja_env = jinja2.Environment(loader=jinja2.FileSystemLoader(template_dir), autoescape=False)

#render a page outside the web framework
def render_string(template,**params):
    t = jinja_env.get_template(template)
    return t.render(params)
def render(template,**kw):
    print(render_string(template,**kw))
#shouldn't be needed anymore but I'll keep it here for future use
#render("welcome_page.html",css=css,js=js)     #redirects to port 8080

#Functions available throughout the framework
class Handler(webapp2.RequestHandler):
    def write(self, *a, **kw):
        self.response.out.write(*a, **kw)      #write plain text

    def render_str(self, template, **params):   #render template as a string
        t = jinja_env.get_template(template)
        return t.render(params)

    def render(self, template, **kw):
        self.write(self.render_str(template, **kw))     #writes the template string

    #Retrieves the running parts_list and session_id if they exist, registers an empty one and an id if they don't
    def get_parts_list(self):
        app = webapp2.get_app()
        session_id = self.request.cookies.get('sessionid')
        if not app.registry.get(session_id):
            app.registry[session_id] = {}
            app.registry[session_id]['parts_list'] = []
        if session_id in app.registry.keys() and session_id is not None:
            parts_list = app.registry.get(session_id)['parts_list']
        else:
            parts_list = []
            session_id = make_session_id()
            self.response.headers.add_header("Set-Cookie","sessionid="+session_id+"; Path=/")
        return parts_list,session_id

    def update_part_list(self,updated_parts_list=None):
        application = webapp2.get_app()
        if updated_parts_list==None:
            exported_map = self.request.POST.get("export_dom")[:-1].split("|")
            if exported_map == ['']:
                parts_list,session_id = self.get_parts_list()
                return parts_list,session_id
            updated_parts_list = []
            for part in exported_map:
                part = tuple(part.split(","))
                if part[1] == "MultiPart":
                    parts_list,session_id = self.get_parts_list()
                    for thingy in parts_list:
                        if thingy.name == part[0]:
                            updated_parts_list.append(thingy)
                else:
                    if part[3]!="":
                        e_part = Part(part[0],part[1],part[2])
                        e_part.description = part[3]
                    else:
                        e_part = Part(part[0],part[1],part[2])
                    updated_parts_list.append(e_part)
        else:
            updated_parts_list = updated_parts_list
        session_id = self.get_parts_list()[1]
        application.registry[session_id]['parts_list'] = updated_parts_list
        return updated_parts_list,session_id

'''
Development note: the webapp2 documentation specifies that retrieving information from the
web page uses the syntax "self.request.get()". While running in its own server, it seems that
different syntax must be used. Use "self.request.POST.get()" for post method data, and
"self.request.GET.get()" for get method data
'''

#Basic object for a plasmid part (description optional)
class Part:
    def __init__(self,name,type,sequence):
        self.name = name
        self.sequence = sequence
        self.type = type
    description = ""
    primer_forward = ""
    primer_reverse = ""
    assembly_method = ""
    bridge_with_next_part = ""
    bridge_with_previous_part = ""

class MultiPart:
    def __init__(self,name,parts):
        self.name = name
        self.parts = parts
        self.type = "MultiPart"
    assembly_method = ""

#Creates a unique session_id for a user to be stored as a cookie. Don't get hungry.
def make_session_id():
    id_nouns = {1:'Cheeseburger',2:'HotDog',3:'CrispyBacon',4:'StrawberryYogurt',5:'SourAppleLollipop',
    6:'HotSauce',7:'MeatballSub',8:'ScrambledEggs',9:'CocoaPuffs',10:'ChocolateMousse',11:'BunchOfBananas',
    12:'FishAndChips',13:'BeefPotRoast',14:'RedBeansAndRice',15:'FriedTofu',16:'EspressoShot',
    17:'CaramelMacchiato',18:'ChocolateChipFrappuccino',19:'HawaiianPizza',20:'RoastedPistachios',
    21:'ChickenTenders',22:'IchirakuRamen',23:'FunnelCake',24:'CottonCandy',25:'SkirtSteak',
    26:'LambChop',27:'ClamChowder',28:'ChickenNoodle',29:'MatzoBall',30:'PotatoChips'
    }
    id_verbs = {1:'Delicious',2:'Appetizing',3:'Delectable',4:'Delightful',5:'Distinctive',
    6:'Enjoyable',7:'Enticing',8:'Exquisite',9:'Heavenly',10:'Luscious',11:'Piquant',
    12:'Pleasant',13:'Rich',14:'Savory',15:'Spicy',16:'Sweet',17:'Tasty',18:'Tempting',
    19:'Yummy',20:'Divine',21:'MouthWatering',22:'Palatable',23:'Sapid',24:'Scrumptious',25:'Tasteful'
                }
    from random import randrange
    import time
    verb = id_verbs[randrange(1,26)]
    noun = id_nouns[randrange(1,31)]
    number = str(randrange(0,999999999999))
    return verb+noun+number

def generateSBOLdoc(builds_list,session_id):
    roles_dict = {"Promoter":"http://identifiers.org/so/SO:0000167",
                  "CDS":"http://identifiers.org/so/SO:0000316",
                  "Terminator":"http://identifiers.org/so/SO:0000141",
                  "RBS":"http://identifiers.org/so/SO:0000139",
                  "oriR":"http://identifiers.org/so/SO:0000296",
                  "userDefined":"http://identifiers.org/so/SO:0000101"
                  }
    import xml.etree.ElementTree as ET
    rdf_root = ET.Element("rdf:RDF",attrib={"xmlns:pr":"http://partsregistry.org","xmlns:rdf":"http://www.w3.org/1999/02/22-rdf-syntax-ns#","xmlns:dcterms":"http://purl.org/dc/terms/","xmlns:prov":"http://www.w3.org/ns/prov#","xmlns:sbol":"http://sbols.org/v2#"})
    xml = ET.ElementTree(element=rdf_root)
    for unpacked_list in builds_list:
        collection_uri = ""
        for i in range(len(unpacked_list)-1):
            collection_uri += unpacked_list[i].name
            collection_uri += "-"
        collection_uri += unpacked_list[-1].name
        collection_node = ET.Element("sbol:Collection",attrib={"rdf:about":"http://ibiocad.igb.illinois.edu/collection/"+collection_uri})
        for part in unpacked_list:
            part_seq_node = ET.Element("sbol:Sequence",attrib={"rdf:about":"http://ibiocad.igb.illinois.edu/seq/"+part.name})
            part_seq_node.append(ET.Element("sbol:persistentIdentity",attrib={"rdf:resource":"http://ibiocad.igb.illinois.edu/seq/"+part.name}))
            seq_disp_id = ET.Element("sbol:displayId",attrib={})
            seq_disp_id.text = part.name
            part_seq_node.append(seq_disp_id)
            sbol_element = ET.Element("sbol:elements",attrib={})
            sbol_element.text = part.sequence
            part_seq_node.append(sbol_element)
            part_seq_node.append(ET.Element("sbol:encoding",attrib={"rdf:resource":"http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html"}))
            collection_node.append(part_seq_node)
            part_node = ET.Element("sbol:ComponentDefinition",attrib={"rdf:about":"http://ibiocad.igb.illinois.edu/component/"+part.name})
            part_node.append(ET.Element("sbol:persistentIdentity",attrib={"rdf:resource":"http://ibiocad.igb.illinois.edu/component/"+part.name}))
            part_disp_id = ET.Element("sbol:displayId",attrib={})
            part_disp_id.text = part.name
            part_node.append(part_disp_id)
            dcterms = ET.Element("dcterms:title",attrib={})
            dcterms.text = part.name
            part_node.append(dcterms)
            dcterms_desc = ET.Element("dcterms:description",attrib={})
            dcterms_desc.text = part.description
            part_node.append(dcterms_desc)
            part_node.append(ET.Element("sbol:type",attrib={"rdf:resource":"http://www.biopax.org/release/biopax-level3.owl#DnaRegion"}))
            part_node.append(ET.Element("sbol:role",attrib={"rdf:resource":roles_dict[part.type]}))
            part_node.append(ET.Element("sbol:sequence",attrib={"rdf:resource":"http://ibiocad.igb.illinois.edu/seq/"+part.name}))
            collection_node.append(part_node)
        rdf_root.append(collection_node)
    xml.write("/var/www/ibiocad/iBioCAD/constructs/plasmid_construct_%s.xml"%session_id)
def reverse_complement(sequence):
    rev_comp = ""
    Watson_Crick = {"A":"T","C":"G","T":"A","G":"C","a":"t","t":"a","c":"g","g":"c"}
    for base in sequence:
        rev_comp = Watson_Crick[base] + rev_comp
    return rev_comp

def create_LCR_bridge(sequence,next_sequence):
    bridge = ""
    current_bridge = ""
    best_tm = -9999999999
    bridge_dict = {}
    for i in range(1,(min(len(sequence),len(next_sequence))+1)):
        current_bridge = reverse_complement(sequence[-i]) + current_bridge
        current_bridge = current_bridge + reverse_complement(next_sequence[i-1])
        bridge_tm = (mt.Tm_NN(current_bridge))
        if bridge_tm == 70.0:
            return current_bridge
        if bridge_tm > 69 and bridge_tm < 71:
            if bridge != "" and abs(70-bridge_tm) < abs(70-best_tm):
                bridge = current_bridge
                best_tm = bridge_tm
        bridge_dict[bridge_tm] = current_bridge
    if bridge == "":
        for tm in bridge_dict.keys():
            if abs(70-tm) < abs(70-best_tm):
                bridge = bridge_dict[tm]
                best_tm = tm
    return bridge

def builds(parts_list):
    import copy
    builds_list = []
    builds = 1
    for part in parts_list:
        if isinstance(part,MultiPart):
            builds *= len(part.parts)
    for i in range(builds):
        build = []
        for part in parts_list:
            if isinstance(part,MultiPart):
                build.append(copy.copy(part.parts[i%len(part.parts)]))
            else:
                build.append(copy.copy(part))
        builds_list.append(build)
    return builds_list

#optimizes the overhangs used in Golden Gate assembly
def golden_gate_optimization(parts_list,backbone_sequence):
    seq_pairs = {0:"ccct",1:"gctc",2:"cggt",3:"gtgc",4:"agcg",5:"ctgt",6:"tgct",7:"atgg",8:"gact",9:"ggac",10:"tccg",11:"ccag",12:"cagc",13:"gttg",14:"cgaa",15:"ccat"}
    seq_matches = []
    for x in range(len(parts_list)+1):
        seq_matches.append([])
        if x == 0:
            for y in range(16):
                if seq_pairs[y] in reverse_complement(backbone_sequence[-35:]):
                    seq_matches[x].append(seq_pairs[y])
                elif seq_pairs[y] in parts_list[x].sequence[:35]:
                    seq_matches[x].append(seq_pairs[y])
        elif x == len(parts_list):
            for y in range(16):
                if seq_pairs[y] in reverse_complement(parts_list[x-1].sequence[-35:]):
                    seq_matches[x].append(seq_pairs[y])
                elif seq_pairs[y] in backbone_sequence[:35]:
                    seq_matches[x].append(seq_pairs[y])
        else:
            for y in range(16):
                if seq_pairs[y] in reverse_complement(parts_list[x].sequence[-35:]):
                    seq_matches[x].append(seq_pairs[y])
                elif seq_pairs[y] in parts_list[x].sequence[:35]:
                    seq_matches[x].append(seq_pairs[y])
    combs = []
    for x in itertools.product(*seq_matches):
        combs.append(x)
    for comb in combs:
        if len(comb) == len(set(comb)):
            return comb
    #if there are no possible combinations
    return None

# '/' path
class MainHandler(Handler):
    def get(self):
        parts_list,session_id = self.get_parts_list()
        self.render("main_page.html",golden_gate_error="",css=css,js=js,parts_list=parts_list)
    def post(self):
        #If "Delete construct" clicked, deletes existing cookie and assigns a new one
        if self.request.POST.get("delete_map")=="TERMINATE":
            self.response.delete_cookie("sessionid")
            self.redirect("/")

        parts_list,session_id = self.get_parts_list()
        if self.request.POST.get("export_dom"):
            parts_list,session_id = self.update_part_list()
        if self.request.POST.get("file_input")!= b'' and self.request.POST.get("file_input") is not None:
            parts_list,session_id = self.get_parts_list()
            file = self.request.POST.get("file_input").file.read().decode("UTF-8")
            if len(list(SeqIO.parse(io.StringIO(file),"fasta"))) > 1:
                application = webapp2.get_app()
                for record in SeqIO.parse(io.StringIO(file),"fasta"):
                    parts_list.append(Part(record.name,"userDefined",record.seq))
                application.registry[session_id]['parts_list'] = parts_list
                self.redirect("/")
            elif len(list(SeqIO.parse(io.StringIO(file),"fasta"))) == 1:
                application = webapp2.get_app()
                for record in SeqIO.parse(io.StringIO(file),"fasta"):
                    name = record.name
                    description = record.name
                    sequence = record.seq
                application.registry[session_id]["file_input"] = [name,sequence,description]
                self.redirect("/inputpart")
            else:
                name = file.split("\r\n")[0][1:]
                description = file.split("\r\n")[0][1:]
                sequence = file.split("\r\n")[1]
                session_id = self.get_parts_list()[1]
                app = webapp2.get_app()
                app.registry[session_id]["file_input"] = [name,sequence,description]
                self.redirect("/inputpart")
        #If "Input New Part" clicked, opens input html page
        if self.request.POST.get("input_part")=="goto_input_window":
            self.redirect("/inputpart")
        if self.request.POST.get("export_map")=="export_map":
            parts_list,session_id = self.update_part_list()
            #for future: data URIs
            #data:[media_type](such as text/plain);charset=utf-8,<data>(encoded in octet characters)
            #sbol uses xml or the python module "sbol" can be used, see docs
        if self.request.POST.get("save_construct") == "save_construct":
            parts_list,session_id = self.update_part_list()
            builds_list = builds(parts_list)
            generateSBOLdoc(builds_list,session_id)
            self.redirect("/construct_download")
        if self.request.POST.get("xml_input")!= b'' and self.request.POST.get("xml_input") is not None:
            rev_roles_dict = {"http://identifiers.org/so/SO:0000167":"Promoter",
                  "http://identifiers.org/so/SO:0000316":"CDS",
                  "http://identifiers.org/so/SO:0000141":"Terminator",
                  "http://identifiers.org/so/SO:0000139":"RBS",
                  "http://identifiers.org/so/SO:0000296":"oriR",
                  "http://identifiers.org/so/SO:0000101":"userDefined"
                  }
            parts_list,session_id = self.get_parts_list()
            xml = self.request.POST.get("xml_input").file.read().decode("UTF-8")
            import xml.etree.ElementTree as ET
            tree = ET.parse(io.StringIO(xml))
            root = tree.getroot()
            builds_list = []
            for unpacked_list_collection in root:
                unpacked_list = []
                components = []
                sequences = []
                for node in unpacked_list_collection:
                    if node.tag == '''{http://sbols.org/v2#}ComponentDefinition''':
                        components.append(node)
                    if node.tag == '''{http://sbols.org/v2#}Sequence''':
                        sequences.append(node)
                for component,sequence in zip(components,sequences):
                    if component[6].get('rdf:resource') == sequence[0].get('rdf:resource'):
                        name = component[1].text
                        type = rev_roles_dict[component[5].get('{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource')]
                        seq = sequence[2].text
                        unpacked_list.append(Part(name,type,seq))
                builds_list.append(unpacked_list)
            if len(builds_list) == 1:
                parts_list,session_id = self.update_part_list(updated_parts_list=builds_list[0])
            else:
                parts_list = []
                for i in range(len(builds_list[0])):
                    is_multipart = False
                    current_part = builds_list[0][i].name
                    for unpacked_list in builds_list:
                        if unpacked_list[i].name != current_part:
                            is_multipart = True
                    if is_multipart:
                        parts = []
                        part_names = []
                        for unpacked_list in builds_list:
                            if unpacked_list[i].name not in part_names:
                                parts.append(unpacked_list[i])
                                part_names.append(unpacked_list[i].name)
                        dynname = ""
                        for i in range(len(parts)-1):
                            dynname += parts[i].name
                            dynname += "--"
                        dynname += parts[-1].name
                        parts_list.append(MultiPart(dynname,parts))
                    else:
                        parts_list.append(builds_list[0][i])
                parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
            self.redirect('/')

        if is_local:
            for record in SeqIO.parse("templates/pET-26b.fa","fasta"):
                default_backbone = record
        else:
            for record in SeqIO.parse("/var/www/ibiocad/iBioCAD/templates/pET-26b.fa","fasta"):
                default_backbone = record
        default_config = {"backbone":default_backbone,"Golden_gate_method":"regular_assembly"}
        parts_list,session_id = self.get_parts_list()
        application = webapp2.get_app()
        if "assembly_config" not in application.registry.get(session_id).keys() or application.registry[session_id]["assembly_config"] is None:
            application.registry[session_id]["assembly_config"] = default_config
        assembly_config = application.registry.get(session_id)["assembly_config"]
        backbone_sequence = assembly_config["backbone"].seq
        new_backbone_sequence = backbone_sequence
        application.registry[session_id]["new_backbone_sequence"] = new_backbone_sequence
        golden_gate_method = assembly_config["Golden_gate_method"]
        #Run Yeast Assembly
        if self.request.POST.get("assembly_method") == "Yeast_Assembly":
            parts_list,session_id = self.update_part_list()
            builds_list = builds(parts_list)
            app = webapp2.get_app()
            app.registry[session_id]['builds_list'] = builds_list
            for unpacked_list in builds_list:
                for part in unpacked_list:
                    part.assembly_method = "Yeast_Assembly"
                for i in range(len(unpacked_list)):
                    if len(unpacked_list) < 2:
                        self.redirect('/')
                    if i == 0:
                        if len(unpacked_list[i].sequence) >= 20:
                            unpacked_list[i].primer_forward = backbone_sequence[-40:] + unpacked_list[i].sequence[:20]
                        else:
                            unpacked_list[i].primer_forward = backbone_sequence[-40:] + unpacked_list[i].sequence
                        if len(unpacked_list[i].sequence) >= 20 and len(unpacked_list[i+1].sequence) >= 40:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-20:] + unpacked_list[i+1].sequence[:40])
                        elif len(unpacked_list[i].sequence) >= 20 and len(unpacked_list[i+1].sequence) < 40:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-20:] + unpacked_list[i+1].sequence)
                        elif len(unpacked_list[i].sequence) < 20 and len(unpacked_list[i+1].sequence) >= 40:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + unpacked_list[i+1].sequence[:40])
                        else:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + unpacked_list[i+1].sequence)
                    elif i == len(unpacked_list)-1:
                        if len(unpacked_list[i].sequence) >= 20:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-20:] + backbone_sequence[:40])
                        else:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + backbone_sequence[:40])
                        if len(unpacked_list[i].sequence) >= 20 and len(unpacked_list[i-1].sequence) >= 40:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence[-40:] + unpacked_list[i].sequence[:20]
                        elif len(unpacked_list[i].sequence) >= 20 and len(unpacked_list[i-1].sequence) < 40:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence + unpacked_list[i].sequence[:20]
                        elif len(unpacked_list[i].sequence) < 20 and len(unpacked_list[i-1].sequence) >= 40:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence[-40:] + unpacked_list[i].sequence
                        else:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence + unpacked_list[i].sequence
                    else:
                        if len(unpacked_list[i].sequence) >= 20 and len(unpacked_list[i-1].sequence) >= 40:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence[-40:] + unpacked_list[i].sequence[:20]
                        elif len(unpacked_list[i].sequence) >= 20 and len(unpacked_list[i-1].sequence) < 40:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence + unpacked_list[i].sequence[:20]
                        elif len(unpacked_list[i].sequence) < 20 and len(unpacked_list[i-1].sequence) >= 40:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence[-40:] + unpacked_list[i].sequence
                        else:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence + unpacked_list[i].sequence
                        if len(unpacked_list[i].sequence) >= 20 and len(unpacked_list[i+1].sequence) >= 40:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-20:] + unpacked_list[i+1].sequence[:40])
                        elif len(unpacked_list[i].sequence) >= 20 and len(unpacked_list[i+1].sequence) < 40:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-20:] + unpacked_list[i+1].sequence)
                        elif len(unpacked_list[i].sequence) < 20 and len(unpacked_list[i+1].sequence) >= 40:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + unpacked_list[i+1].sequence[:40])
                        else:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + unpacked_list[i+1].sequence)
                #parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
                for i in range(len(unpacked_list)):
                    if len(unpacked_list) < 2:
                        self.redirect('/')
                    if i == 0:
                        if len(unpacked_list[i+1].sequence) >= 40:
                            unpacked_list[i].sequence = backbone_sequence[-40:] + unpacked_list[i].sequence + unpacked_list[i+1].sequence[0:40]
                        else:
                            unpacked_list[i].sequence = backbone_sequence[-40:] + unpacked_list[i].sequence + unpacked_list[i+1].sequence
                    elif i == len(unpacked_list)-1:
                        if len(unpacked_list[i-1].sequence) >= 40:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence[-80:-40] + unpacked_list[i].sequence + backbone_sequence[:40]
                        else:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence + unpacked_list[i].sequence + backbone_sequence[:40]
                    else:
                        if len(unpacked_list[i-1].sequence) >= 40 and len(unpacked_list[i+1].sequence) >= 40:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence[-80:-40] + unpacked_list[i].sequence + unpacked_list[i+1].sequence[0:40]
                        elif len(unpacked_list[i-1].sequence) >= 40:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence[-80:-40] + unpacked_list[i].sequence + unpacked_list[i+1].sequence
                        elif len(unpacked_list[i+1].sequence) >= 40:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence + unpacked_list[i].sequence + unpacked_list[i+1].sequence[0:40]
                        else:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence + unpacked_list[i].sequence + unpacked_list[i+1].sequence
            parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
            app.registry[session_id]['builds_list'] = builds_list
            self.redirect("/assembly")
        #Run Gibson Assembly
        if self.request.POST.get("assembly_method") == "Gibson_Assembly":
            parts_list,session_id = self.update_part_list()
            builds_list = builds(parts_list)
            app = webapp2.get_app()
            app.registry[session_id]['builds_list'] = builds_list
            for unpacked_list in builds_list:
                for part in unpacked_list:
                    part.assembly_method = "Gibson_Assembly"
                for i in range(len(unpacked_list)):
                    if len(unpacked_list) < 2:
                        self.redirect('/')
                    if i == 0:
                        if len(unpacked_list[i].sequence) >= 25:
                            unpacked_list[i].primer_forward = backbone_sequence[-25:] + unpacked_list[i].sequence[:25]
                        else:
                            unpacked_list[i].primer_forward = backbone_sequence[-25:] + unpacked_list[i].sequence
                        if len(unpacked_list[i].sequence) >= 25 and len(unpacked_list[i+1].sequence) >= 25:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-25:] + unpacked_list[i+1].sequence[:25])
                        elif len(unpacked_list[i].sequence) >= 25 and len(unpacked_list[i+1].sequence) < 25:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-25:] + unpacked_list[i+1].sequence)
                        elif len(unpacked_list[i].sequence) < 25 and len(unpacked_list[i+1].sequence) >= 25:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + unpacked_list[i+1].sequence[:25])
                        else:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + unpacked_list[i+1].sequence)
                    elif i == len(unpacked_list)-1:
                        if len(unpacked_list[i].sequence) >= 25:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-25:] + backbone_sequence[:25])
                        else:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + backbone_sequence[:25])
                        if len(unpacked_list[i].sequence) >= 25 and len(unpacked_list[i-1].sequence) >= 25:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence[-25:] + unpacked_list[i].sequence[:25]
                        elif len(unpacked_list[i].sequence) >= 25 and len(unpacked_list[i-1].sequence) < 25:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence + unpacked_list[i].sequence[:25]
                        elif len(unpacked_list[i].sequence) < 25 and len(unpacked_list[i-1].sequence) >= 25:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence[-25:] + unpacked_list[i].sequence
                        else:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence + unpacked_list[i].sequence
                    else:
                        if len(unpacked_list[i].sequence) >= 25 and len(unpacked_list[i-1].sequence) >= 25:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence[-25:] + unpacked_list[i].sequence[:25]
                        elif len(unpacked_list[i].sequence) >= 25 and len(unpacked_list[i-1].sequence) < 25:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence + unpacked_list[i].sequence[:25]
                        elif len(unpacked_list[i].sequence) < 25 and len(unpacked_list[i-1].sequence) >= 25:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence[-25:] + unpacked_list[i].sequence
                        else:
                            unpacked_list[i].primer_forward = unpacked_list[i-1].sequence + unpacked_list[i].sequence
                        if len(unpacked_list[i].sequence) >= 25 and len(unpacked_list[i+1].sequence) >= 25:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-25:] + unpacked_list[i+1].sequence[:25])
                        elif len(unpacked_list[i].sequence) >= 25 and len(unpacked_list[i+1].sequence) < 25:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-25:] + unpacked_list[i+1].sequence)
                        elif len(unpacked_list[i].sequence) < 25 and len(unpacked_list[i+1].sequence) >= 25:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + unpacked_list[i+1].sequence[:25])
                        else:
                            unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence + unpacked_list[i+1].sequence)
                #parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
                for i in range(len(unpacked_list)):
                    if len(unpacked_list) < 2:
                        self.redirect('/')
                    if i == 0:
                        if len(unpacked_list[i+1].sequence) >= 25:
                            unpacked_list[i].sequence = backbone_sequence[-25:] + unpacked_list[i].sequence + unpacked_list[i+1].sequence[0:25]
                        else:
                            unpacked_list[i].sequence = backbone_sequence[-25:] + unpacked_list[i].sequence + unpacked_list[i+1].sequence
                    elif i == len(unpacked_list)-1:
                        if len(unpacked_list[i-1].sequence) >= 25:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence[-50:-25] + unpacked_list[i].sequence + backbone_sequence[:25]
                        else:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence + unpacked_list[i].sequence + backbone_sequence[:25]
                    else:
                        if len(unpacked_list[i-1].sequence) >= 25 and len(unpacked_list[i+1].sequence) >= 25:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence[-50:-25] + unpacked_list[i].sequence + unpacked_list[i+1].sequence[0:25]
                        elif len(unpacked_list[i-1].sequence) >= 25:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence[-50:-25] + unpacked_list[i].sequence + unpacked_list[i+1].sequence
                        elif len(unpacked_list[i+1].sequence) >= 25:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence + unpacked_list[i].sequence + unpacked_list[i+1].sequence[0:25]
                        else:
                            unpacked_list[i].sequence = unpacked_list[i-1].sequence + unpacked_list[i].sequence + unpacked_list[i+1].sequence
            parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
            app.registry[session_id]['builds_list'] = builds_list
            self.redirect("/assembly")
        if self.request.POST.get("assembly_method") == "LCR":
            parts_list,session_id = self.update_part_list()
            builds_list = builds(parts_list)
            application = webapp2.get_app()
            application.registry[session_id]['builds_list'] = builds_list
            for unpacked_list in builds_list:
                for part in unpacked_list:
                    part.assembly_method = "LCR"
                for i in range(len(unpacked_list)):
                    if len(unpacked_list)<2:
                        self.redirect('/')
                    if i == (len(unpacked_list)-1):
                        unpacked_list[i].bridge_with_next_part = create_LCR_bridge(unpacked_list[i].sequence,backbone_sequence[:200])
                    if i == 0:
                        unpacked_list[i].bridge_with_previous_part = create_LCR_bridge(backbone_sequence[-200:],unpacked_list[i].sequence)
                    if i < (len(unpacked_list)-1):
                        unpacked_list[i].bridge_with_next_part = create_LCR_bridge(unpacked_list[i].sequence,unpacked_list[i+1].sequence)
            parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
            application.registry[session_id]['builds_list'] = builds_list
            self.redirect("/assembly")
        golden_gate_error=""
        if self.request.POST.get("assembly_method") == "Type_II_Restriction_Enzyme":
            parts_list,session_id = self.update_part_list()
            builds_list = builds(parts_list)
            app = webapp2.get_app()
            app.registry[session_id]['builds_list'] = builds_list
            new_backbone_sequence = backbone_sequence
            golden_gate_overhangs = [
                "ccct","gctc","cggt","gtgc","agcg","ctgt","tgct","atgg","gact","ggac","tccg","ccag","cagc","gttg","cgaa","ccat"
            ]
            if golden_gate_method == "regular_assembly":
                for unpacked_list in builds_list:
                    for part in unpacked_list:
                        part.assembly_method = "Type_II_Restriction_Enzyme"
                        if "ggtctc" in part.sequence.lower() or "gagacc" in part.sequence.lower():
                            golden_gate_error = "BsaI_in_seq" +"|"+ part.name + "|" + part.sequence.lower()
                    for i in range(len(unpacked_list)):
                        if len(unpacked_list)<2:
                            self.redirect('/')
                        if golden_gate_error != "":
                            break
                        if len(unpacked_list) > 16:
                            break
                        unpacked_list[i].primer_forward = "aaggtctca" + golden_gate_overhangs[i] + unpacked_list[i].sequence[:20]
                        unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-20:] + golden_gate_overhangs[i+1] + "agagaccaa")
                        unpacked_list[i].sequence = "aaggtctca" + golden_gate_overhangs[i] + unpacked_list[i].sequence + golden_gate_overhangs[i+1] + "agagaccaa"
            if golden_gate_method == "scarless_assembly":
                for unpacked_list in builds_list:
                    for part in unpacked_list:
                        part.assembly_method = "Type_II_Restriction_Enzyme"
                        if "ggtctc" in part.sequence.lower() or "gagacc" in part.sequence.lower():
                            golden_gate_error = "BsaI_in_seq" +"|"+ part.name + "|" + part.sequence.lower()
                    for i in range(len(unpacked_list)):
                        if len(unpacked_list)<2:
                            self.redirect('/')
                        if golden_gate_error != "":
                            break
                        if len(unpacked_list) > 16:
                            break
                        gg_opt = golden_gate_optimization(unpacked_list,backbone_sequence)
                        if gg_opt is None:
                            golden_gate_error = "no_efficient_overhang_combinations"
                            break
                        if i == 0:
                            if gg_opt[0] in new_backbone_sequence[-35:]:
                                new_backbone_sequence = new_backbone_sequence[:new_backbone_sequence.find(gg_opt[0])] + gg_opt[0] + "agagaccaa"
                                unpacked_list[0].primer_forward = "aaggtctca" + new_backbone_sequence[new_backbone_sequence.find(gg_opt[0]):-9] + unpacked_list[0].sequence[:20]
                                unpacked_list[0].sequence = "aaggtctca" + new_backbone_sequence[new_backbone_sequence.find(gg_opt[0]):-9] + unpacked_list[0].sequence
                            elif gg_opt[0] in unpacked_list[0].sequence[:35]:
                                new_backbone_sequence = new_backbone_sequence + unpacked_list[0].sequence[:unpacked_list[0].sequence.find(gg_opt[0])] + gg_opt[0] + "agagaccaa"
                                unpacked_list[0].primer_forward = "aaggtctca" + unpacked_list[0].sequence[unpacked_list[0].sequence.find(gg_opt[0]):unpacked_list[0].sequence.find(gg_opt[0])+20]
                                unpacked_list[0].sequence = "aaggtctca" + unpacked_list[0].sequence[unpacked_list[0].sequence.find(gg_opt[0]):]
                        elif i == len(unpacked_list)-1:
                            if gg_opt[-1] in new_backbone_sequence[:35]:
                                new_backbone_sequence = "aaggtctca" + new_backbone_sequence[new_backbone_sequence.find(gg_opt[-1]):]
                                unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[-1].sequence[-20:] + new_backbone_sequence[9:new_backbone_sequence.find(gg_opt[-1])] + gg_opt[-1] + "agagaccaa")
                                unpacked_list[-1].sequence = unpacked_list[-1].sequence + new_backbone_sequence[:new_backbone_sequence.find(gg_opt[-1])] + gg_opt[-1] + "agagaccaa"
                            elif gg_opt[-1] in unpacked_list[-1].sequence[-35:]:
                                new_backbone_sequence = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-1]):] + new_backbone_sequence
                                unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-1])-20:unpacked_list[-1].sequence.find(gg_opt[-1])] + gg_opt[-1] + "agagaccaa")
                                unpacked_list[-1].sequence = unpacked_list[-1].sequence[:unpacked_list[-1].sequence.find(gg_opt[-1])] + gg_opt[-1] + "agagaccaa"
                            if gg_opt[-2] in unpacked_list[-1].sequence[:35]:
                                unpacked_list[-2].primer_reverse = reverse_complement(unpacked_list[-2].sequence[-20:] + unpacked_list[-1].sequence[:unpacked_list[-1].sequence.find(gg_opt[-2])] + gg_opt[-2] + "agagaccaa")
                                unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-2]):unpacked_list[-1].sequence.find(gg_opt[-2])+20]
                                unpacked_list[-1].sequence = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-2]):]
                                unpacked_list[-2].sequence = unpacked_list[-2].sequence + unpacked_list[-1].sequence[:unpacked_list[-1].sequence.find(gg_opt[-2])] + gg_opt[-2] + "agagaccaa"
                            elif gg_opt[-2] in unpacked_list[-2].sequence[-35:]:
                                unpacked_list[-2].primer_reverse = reverse_complement(unpacked_list[-2].sequence[unpacked_list[-2].sequence.find(gg_opt[-2])-20:unpacked_list[-2].sequence.find(gg_opt[-2])] + gg_opt[-2] + "agagaccaa")
                                unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[-2].sequence[unpacked_list[-2].sequence.find(gg_opt[-2]):] + unpacked_list[-1].sequence[:20]
                                unpacked_list[-1].sequence = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-2]):]
                                unpacked_list[-2].sequence = unpacked_list[-2].sequence[:unpacked_list[-2].sequence.find(gg_opt[-2])] + gg_opt[-2] + "agagaccaa"
                        else:
                            if gg_opt[i] in unpacked_list[i].sequence[:35]:
                                unpacked_list[i-1].primer_reverse = reverse_complement(unpacked_list[i-1].sequence[-20:] + unpacked_list[i].sequence[:unpacked_list[i].sequence.find(gg_opt[i])] + gg_opt[i] + "agagaccaa")
                                unpacked_list[i].primer_forward = "aaggtctca" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(gg_opt[i]):unpacked_list[i].sequence.find(gg_opt[i])+20]
                                unpacked_list[i].sequence = "aaggtctca" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(gg_opt[i]):]
                                unpacked_list[i-1].sequence = unpacked_list[i-1].sequence + unpacked_list[i].sequence[:unpacked_list[i].sequence.find(gg_opt[i])] + gg_opt[i] + "agagaccaa"
                            elif gg_opt[i] in unpacked_list[i-1].sequence[-35:]:
                                unpacked_list[i-1].primer_reverse = reverse_complement(unpacked_list[i-1].sequence[unpacked_list[i-1].sequence.find(gg_opt[i])-20:unpacked_list[i-1].sequence.find(gg_opt[i])] + gg_opt[i] + "agagaccaa")
                                unpacked_list[i].primer_forward = "aaggtctca" + unpacked_list[i-1].sequence[unpacked_list[i-1].sequence.find(gg_opt[i]):] + unpacked_list[i].sequence[:20]
                                unpacked_list[i].sequence = "aaggtctca" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(gg_opt[i]):]
                                unpacked_list[i-1].sequence = unpacked_list[i-1].sequence[:unpacked_list[i-1].sequence.find(gg_opt[i])] + gg_opt[i] + "agagaccaa"
            parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
            app.registry[session_id]['builds_list'] = builds_list
            app.registry[session_id]["new_backbone_sequence"] = new_backbone_sequence
            if golden_gate_error == "":
                self.redirect("/assembly")
        self.render("main_page.html",golden_gate_error=golden_gate_error,css=css,js=js,parts_list=parts_list)


'''
Go to the path "/break" to break the process so that it can (or must) be restarted.
Set server_debug to False to disable this (or it can be deleted entirely).
Break password is "illini".
'''
class breakHandler(Handler):
    def get(self):
        self.render("break_page.html",breakpass_in="",break_error="")
    def post(self):
        #If enabled and correct password is entered, shuts off current process
        breakpass = self.request.POST.get("breakpass")
        if breakpass == "illini" and server_debug:
            pid = os.getpid()
            import signal
            os.kill(pid,signal.SIGINT)
        elif breakpass == "illini":
            self.render("break_page.html",breakpass_in=breakpass,break_error="debug mode is off")
        else:
            self.render("break_page.html",breakpass_in=breakpass,break_error="wrong password")


'''
 '/inputpart'
Gets data from the input page. If the required fields are entered, creates an instance
of the Part class and adds it to the running parts_list if it exists or creates a new
one if not.
'''
class InputPartHandler(Handler):
    def get(self):
        name = ""
        sequence = ""
        description = ""
        session_id = self.get_parts_list()[1]
        application = webapp2.get_app()
        if "file_input" in application.registry.get(session_id).keys():
            file = application.registry.get(session_id).pop("file_input")
            name += file[0]
            sequence += file[1]
            description += file[2]
        self.render("input_part.html",css=css,js=js,nameerror="",sequenceerror="",name=name,sequence=sequence,description=description,dynname="")
    def post(self):
        if self.request.POST.get('cancel')=="cancel":
            self.redirect("/")
        if self.request.POST.get('add_part')=="add":
            if self.request.POST.getall("files")!= b'' and self.request.POST.get("files") is not None and self.request.POST.get('input_type') == 'dynamic':
                dynname = self.request.POST.get('dynname')
                parts = []
                nameerror = ""
                dynnameerror = ""
                fileerror = ""
                files = self.request.POST.getall("files")
                parts_list,session_id = self.get_parts_list()
                for part in parts_list:
                    if dynname == part.name:
                        dynnameerror = "That name has already been input!"
                if dynname == "":
                    dynnameerror = "You must specify a name!"
                if files == [b'']:
                    fileerror = "You must upload a file!"
                if dynnameerror == "" and fileerror == "":
                    for f in files:
                        if len(files) == 0:
                            self.redirect("/")
                        data = f.file.read().decode("UTF-8").split("\r\n")
                        mpart = Part(data[0][1:],"userDefined",data[1])
                        mpart.is_sub_node = True
                        mpart.parentname = dynname
                        parts.append(mpart)
                    multipart = MultiPart(dynname,parts)
                    parts_list,session_id = self.get_parts_list()
                    parts_list.append(multipart)
                    application.registry[session_id]['parts_list'] = parts_list
                    self.redirect("/")
                self.render("input_part.html",css=css,js=js,nameerror="",dynnameerror=dynnameerror,fileerror=fileerror,dynname=dynname,dynamic="yes")

            if self.request.POST.get('input_type') == "library":
                library_imports = self.request.POST.getall("library_inputs")
                parts_list,session_id = self.get_parts_list()
                for library_input in library_imports:
                    for record in SeqIO.parse("/var/www/ibiocad/iBioCAD/templates/%s.txt"%library_input.split("|")[0],"fasta"):
                        parts_list.append(Part(record.name,library_input.split("|")[1],record.seq))
                application.registry[session_id]['parts_list'] = parts_list
                self.redirect("/")

            if self.request.POST.get('input_type') == "static":
                dynname = ""
                name = self.request.POST.get('name')
                type = self.request.POST.get('type')
                sequence = self.request.POST.get('sequence')
                description = self.request.POST.get('description')
                nameerror=""
                sequenceerror=""
                if name == "" or sequence == "":
                    if not name:
                        nameerror = "You must specify a name!"
                    if not sequence:
                        sequenceerror = "You must specify a sequence!"
                parts_list,session_id = self.get_parts_list()
                for part in parts_list:
                    if name == part.name:
                        nameerror = "That name has already been input!"
                for b in sequence:
                    if b not in "AGCTacgt":
                        sequenceerror = "That's not a valid DNA sequence!"
                if nameerror == "" and sequenceerror == "":
                    part = Part(name,type,sequence)
                    if description:
                        part.description = description
                    parts_list,session_id = self.get_parts_list()
                    parts_list.append(part)
                    application.registry[session_id]['parts_list'] = parts_list
                    self.redirect("/")
                self.render("input_part.html",css=css,js=js,nameerror=nameerror,sequenceerror=sequenceerror,name=name,sequence=sequence,type=type,description=description,dynname=dynname)

class AssemblyHandler(Handler):
    def get(self):
        parts_list,session_id = self.get_parts_list()
        application = webapp2.get_app()
        builds_list = application.registry.get(session_id)['builds_list']
        new_backbone_sequence = application.registry.get(session_id)["new_backbone_sequence"]
        import csv
        if is_local:
            assembly_file = "constructs/plasmid_assembly_%s.csv"
        else:
            assembly_file = "/var/www/ibiocad/iBioCAD/constructs/plasmid_assembly_%s.csv"
        with open(assembly_file%session_id,'w') as csvfile:
            fieldnames = ['Name','Type','Sequence','Description']
            csvdictwriter = csv.DictWriter(csvfile,fieldnames=fieldnames)
            csvwriter = csv.writer(csvfile)
            for unpacked_list in builds_list:
                if len(unpacked_list)<2:
                    break
                if unpacked_list[0].assembly_method == "Yeast_Assembly" or unpacked_list[0].assembly_method == "Gibson_Assembly" or unpacked_list[0].assembly_method == "Type_II_Restriction_Enzyme":
                    csvdictwriter2 = csv.DictWriter(csvfile,fieldnames=['Name','','Sequence'])
                elif unpacked_list[0].assembly_method == "LCR":
                    csvdictwriter2 = csv.DictWriter(csvfile,fieldnames=['Bridge','','Sequence'])
                build = []
                for part in unpacked_list:
                    build.append(part.name)
                csvwriter.writerow([','.join(build)])
                csvwriter.writerow(['Modified Parts'])
                try:
                    csvdictwriter.writeheader()
                except:
                    csvwriter.writerow(fieldnames)
                csvwriter.writerow(['Backbone','',new_backbone_sequence,''])
                for part in unpacked_list:
                    csvdictwriter.writerow({'Name':part.name,'Type':part.type,'Sequence':part.sequence,'Description':part.description})
                csvwriter.writerow([])
                if unpacked_list[0].assembly_method == "Yeast_Assembly" or unpacked_list[0].assembly_method == "Gibson_Assembly" or unpacked_list[0].assembly_method == "Type_II_Restriction_Enzyme":
                    csvwriter.writerow(['Primers'])
                elif unpacked_list[0].assembly_method == "LCR":
                    csvwriter.writerow(['Bridges'])
                try:
                    csvdictwriter2.writeheader()
                except:
                    if unpacked_list[0].assembly_method == "Yeast_Assembly" or unpacked_list[0].assembly_method == "Gibson_Assembly" or unpacked_list[0].assembly_method == "Type_II_Restriction_Enzyme":
                        csvwriter.writerow(['Name','','Sequence'])
                    elif unpacked_list[0].assembly_method == "LCR":
                        csvwriter.writerow(['Bridge','','Sequence'])
                if unpacked_list[0].assembly_method == "Yeast_Assembly" or unpacked_list[0].assembly_method == "Gibson_Assembly" or unpacked_list[0].assembly_method == "Type_II_Restriction_Enzyme":
                    for part in unpacked_list:
                        csvdictwriter2.writerow({'Name':(part.name+",forward"),'':'','Sequence':part.primer_forward})
                        csvdictwriter2.writerow({'Name':(part.name+",reverse"),'':'','Sequence':part.primer_reverse})
                elif unpacked_list[0].assembly_method == "LCR":
                    csvdictwriter2.writerow({'Bridge':("plasmid backbone-"+parts_list[0].name),'':'','Sequence':parts_list[0].bridge_with_previous_part})
                    for i,part in enumerate(unpacked_list[:-1]):
                        csvdictwriter2.writerow({'Bridge':(part.name+"-"+unpacked_list[i+1].name),'':'','Sequence':part.bridge_with_next_part})
                    csvdictwriter2.writerow({'Bridge':(parts_list[-1].name+"-plasmid backbone"),'':'','Sequence':parts_list[-1].bridge_with_next_part})
                csvwriter.writerow([])
        with open(assembly_file%session_id,'r') as f:
            data_uri = "data:text/csv;base64,"
            import base64
            if sys.version_info[0] == 3:
                data_uri += str(base64.b64encode(bytes(f.read(),'utf-8')),'utf-8')
            elif sys.version_info[0] == 2:
                data_uri += base64.b64encode(f.read())
        import datetime
        filename = "plasmid_assembly_%s.csv"%datetime.date.today()
        os.remove(assembly_file%session_id)
        part_dir = builds_list[0][0].assembly_method
        self.render("assembly_page.html",builds_list=builds_list,part_dir=part_dir,new_backbone_sequence=new_backbone_sequence,data_uri=data_uri,filename=filename,css=css)
class ConstructDownloadHandler(Handler):
    def get(self):
        parts_list,session_id = self.get_parts_list()
        if is_local:
            construct_file = "constructs/plasmid_construct_%s.xml"
        else:
            construct_file = "/var/www/ibiocad/iBioCAD/constructs/plasmid_construct_%s.xml"
        with open(construct_file%session_id,'r') as construct:
            file = construct.read()
            data_uri = "data:text/xml;base64,"
            import base64
            try:
                data_uri += base64.b64encode(file)
            except:
                data_uri += str(base64.b64encode(bytes(file,'utf-8')),'utf-8')
        import datetime
        filename = "plasmid_construct_%s.xml"%datetime.date.today()
        os.remove(construct_file%session_id)
        self.render("construct_download.html",css=css,data_uri=data_uri,filename=filename)

class ConfigHandler(Handler):
    def get(self):
        if is_local:
            #path on local
            for record in SeqIO.parse("templates/pET-26b.fa","fasta"):
                default_backbone = record
        else:
            #path on server
            for record in SeqIO.parse("/var/www/ibiocad/iBioCAD/templates/pET-26b.fa","fasta"):
                default_backbone = record
        default_config = {"backbone":default_backbone,"Golden_gate_method":"regular_assembly"}
        parts_list,session_id = self.get_parts_list()
        application = webapp2.get_app()
        if "assembly_config" not in application.registry.get(session_id).keys() or application.registry[session_id]["assembly_config"] is None:
            application.registry[session_id]["assembly_config"] = default_config
        assembly_config = application.registry.get(session_id)["assembly_config"]
        self.render("config.html",css=css,assembly_config=assembly_config)
    def post(self):
        parts_list,session_id = self.get_parts_list()
        application = webapp2.get_app()
        if self.request.POST.get('cancel')=="cancel":
            self.redirect("/")
        if self.request.POST.get('save')=="save":
            for key in application.registry.get(session_id)["assembly_config"].keys():
                if key == "backbone":
                    if self.request.POST.get("backbone")!= b'' and self.request.POST.get("backbone") is not None:
                        file = self.request.POST.get("backbone").file.read().decode("UTF-8")#.split("\r\n")
                        import io
                        textfile = io.StringIO(file)
                        for record in SeqIO.parse(textfile,"fasta"):
                            application.registry[session_id]["assembly_config"][key] = record
                        textfile.close()
                    else:
                        if is_local:
                            #path on local
                            for record in SeqIO.parse("templates/pET-26b.fa","fasta"):
                                application.registry[session_id]["assembly_config"][key] = record
                        else:
                            #path on server
                            for record in SeqIO.parse("/var/www/ibiocad/iBioCAD/templates/pET-26b.fa","fasta"):
                                application.registry[session_id]["assembly_config"][key] = record
                else:
                    application.registry[session_id]["assembly_config"][key] = self.request.POST.get(key)
            self.redirect("/")

#Initialize the web framework. URL paths are attached to their respective handlers
#WSGI app must have the variable name "application" I think
application = webapp2.WSGIApplication([
    ('/', MainHandler),('/break',breakHandler),('/inputpart',InputPartHandler),('/assembly',AssemblyHandler),('/construct_download',ConstructDownloadHandler),('/config',ConfigHandler)#,('/ExamplePath',exampleHandler)
], debug=True)
application.set_globals(app=application)

#Serve the web app on a specific IP and port, useful for local deployment
#Additional helpful arguments available in httpserver documentation
def main():
    from paste import httpserver
    if is_local:
        httpserver.serve(app, host="0.0.0.0",port="80")
    else:
        pass

if __name__ == "__main__":
    main()