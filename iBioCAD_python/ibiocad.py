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

def generateSBOLdoc(parts_list,session_id):
    roles_dict = {"Promoter":"http://identifiers.org/so/SO:0000167",
                  "CDS":"http://identifiers.org/so/SO:0000316",
                  "Terminator":"http://identifiers.org/so/SO:0000141",
                  "RBS":"http://identifiers.org/so/SO:0000139",
                  "oriR":"http://identifiers.org/so/SO:0000296",
                  "userDefined":"http://identifiers.org/so/SO:0000101"
                  }
    xml = """<?xml version="1.0" ?>
    <rdf:RDF xmlns:pr="http://partsregistry.org" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns:dcterms="http://purl.org/dc/terms/" xmlns:prov="http://www.w3.org/ns/prov#" xmlns:sbol="http://sbols.org/v2#">\n"""
    for part in parts_list:
        if isinstance(part,MultiPart):
            for mpart in part.parts:
                xml += """<sbol:Sequence rdf:about="http://ibiocadsite.web.engr.illinois.edu/seq/""" + mpart.name + '">\n'
                xml += """\t<sbol:persistentIdentity rdf:resource="http://ibiocadsite.web.engr.illinois.edu/seq/""" + mpart.name + '"/>\n'
                xml += """\t<sbol:displayId>""" + mpart.name + """</sbol:displayId>\n"""
                xml += """\t<sbol:elements>""" + mpart.sequence + """</sbol:elements>\n"""
                xml += """\t<sbol:encoding rdf:resource="http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html"/>\n"""
                xml += """</sbol:Sequence>\n"""
                xml += """<sbol:ComponentDefinition rdf:about="http://ibiocadsite.web.engr.illinois.edu/""" + mpart.name + '">\n'
                xml += """\t<sbol:persistentIdentity rdf:resource="http://ibiocadsite.web.engr.illinois.edu/""" + mpart.name + '"/>\n'
                xml += """\t<sbol:displayId>""" + mpart.name + """</sbol:displayId>\n"""
                xml += """\t<dcterms:title>""" + mpart.name + """</dcterms:title>\n"""
                xml += """\t<dcterms:description>""" + mpart.description + """</dcterms:description>\n"""
                xml += """\t<sbol:type rdf:resource="http://www.biopax.org/release/biopax-level3.owl#DnaRegion"/>\n"""
                xml += '''\t<sbol:role rdf:resource="''' + roles_dict[part.type] + '"/>\n'
                xml += """\t<sbol:sequence rdf:resource="http://ibiocadsite.web.engr.illinois.edu/seq/""" + mpart.name + '"/>\n'
                xml += """</sbol:ComponentDefinition>\n"""
        xml += """<sbol:Sequence rdf:about="http://ibiocadsite.web.engr.illinois.edu/seq/""" + part.name + '">\n'
        xml += """\t<sbol:persistentIdentity rdf:resource="http://ibiocadsite.web.engr.illinois.edu/seq/""" + part.name + '"/>\n'
        xml += """\t<sbol:displayId>""" + part.name + """</sbol:displayId>\n"""
        xml += """\t<sbol:elements>""" + part.sequence + """</sbol:elements>\n"""
        xml += """\t<sbol:encoding rdf:resource="http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html"/>\n"""
        xml += """</sbol:Sequence>\n"""
        xml += """<sbol:ComponentDefinition rdf:about="http://ibiocadsite.web.engr.illinois.edu/""" + part.name + '">\n'
        xml += """\t<sbol:persistentIdentity rdf:resource="http://ibiocadsite.web.engr.illinois.edu/""" + part.name + '"/>\n'
        xml += """\t<sbol:displayId>""" + part.name + """</sbol:displayId>\n"""
        xml += """\t<dcterms:title>""" + part.name + """</dcterms:title>\n"""
        xml += """\t<dcterms:description>""" + part.description + """</dcterms:description>\n"""
        xml += """\t<sbol:type rdf:resource="http://www.biopax.org/release/biopax-level3.owl#DnaRegion"/>\n"""
        xml += '''\t<sbol:role rdf:resource="''' + roles_dict[part.type] + '"/>\n'
        xml += """\t<sbol:sequence rdf:resource="http://ibiocadsite.web.engr.illinois.edu/seq/""" + part.name + '"/>\n'
        xml += """</sbol:ComponentDefinition>\n"""
    xml += """</rdf:RDF>"""
    with open("constructs/plasmid_construct_%s.xml"%session_id,"w") as construct:
        construct.write(xml)

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
            file = self.request.POST.get("file_input").file.read().decode("UTF-8").split("\r\n")
            name = file[0][1:]
            description = file[0][1:]
            sequence = file[1]
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
            generateSBOLdoc(parts_list,session_id)
            self.redirect("/construct_download")
        if is_local:
            for record in SeqIO.parse("templates/pET-26b.fa","fasta"):
                default_backbone = record
        else:
            for record in SeqIO.parse("/var/www/ibiocad/iBioCAD/templates/pET-26b.fa","fasta"):
                default_backbone = record
        default_config = {"backbone":default_backbone}
        parts_list,session_id = self.get_parts_list()
        app = webapp2.get_app()
        if "assembly_config" not in app.registry.get(session_id).keys() or app.registry[session_id]["assembly_config"] is None:
            app.registry[session_id]["assembly_config"] = default_config
        assembly_config = app.registry.get(session_id)["assembly_config"]
        backbone_sequence = assembly_config["backbone"].seq
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
                        break
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
                        break
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
                        break
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
                        break
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
            self.redirect("/assembly")
        if self.request.POST.get("assembly_method") == "LCR":
            parts_list,session_id = self.update_part_list()
            builds_list = builds(parts_list)
            app = webapp2.get_app()
            app.registry[session_id]['builds_list'] = builds_list
            for unpacked_list in builds_list:
                for part in unpacked_list:
                    part.assembly_method = "LCR"
                for i in range(len(unpacked_list)):
                    if len(unpacked_list)<2:
                        break
                    if i == (len(unpacked_list)-1):
                        unpacked_list[i].bridge_with_next_part = create_LCR_bridge(unpacked_list[i].sequence,backbone_sequence[:200])
                    if i == 0:
                        unpacked_list[i].bridge_with_previous_part = create_LCR_bridge(backbone_sequence[-200:],unpacked_list[i].sequence)
                    if i < (len(unpacked_list)-1):
                        unpacked_list[i].bridge_with_next_part = create_LCR_bridge(unpacked_list[i].sequence,unpacked_list[i+1].sequence)
            parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
            self.redirect("/assembly")
        golden_gate_error=""
        if self.request.POST.get("assembly_method") == "Type_II_Restriction_Enzyme":
            parts_list,session_id = self.update_part_list()
            builds_list = builds(parts_list)
            app = webapp2.get_app()
            app.registry[session_id]['builds_list'] = builds_list
            golden_gate_overhangs = [
                "ccct","gctc","cggt","gtgc","agcg","ctgt","tgct","atgg","gact","ggac","tccg","ccag","cagc","gttg","cgaa","ccat"
            ]
            for unpacked_list in builds_list:
                for part in unpacked_list:
                    part.assembly_method = "Type_II_Restriction_Enzyme"
                    if "ggtctc" or "gagacc" in part.sequence.lower():
                        golden_gate_error = "BsaI_in_seq"
                for i in range(len(unpacked_list)):
                    if golden_gate_error != "":
                        break
                    if len(unpacked_list) > 16:
                        break
                    unpacked_list[i].primer_forward = "aaggtctca" + golden_gate_overhangs[i] + unpacked_list[i].sequence[:20]
                    unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-20:] + golden_gate_overhangs[i+1] + "agagaccaa")
                    unpacked_list[i].sequence = "aaggtctca" + golden_gate_overhangs[i] + unpacked_list[i].sequence + golden_gate_overhangs[i+1] + "agagaccaa"
            parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
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
        app = webapp2.get_app()
        if "file_input" in app.registry.get(session_id).keys():
            file = app.registry.get(session_id).pop("file_input")
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
        app = webapp2.get_app()
        builds_list = app.registry.get(session_id)['builds_list']
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
        self.render("assembly_page.html",builds_list=builds_list,data_uri=data_uri,filename=filename,css=css)

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
        default_config = {"backbone":default_backbone}
        parts_list,session_id = self.get_parts_list()
        app = webapp2.get_app()
        if "assembly_config" not in app.registry.get(session_id).keys() or app.registry[session_id]["assembly_config"] is None:
            app.registry[session_id]["assembly_config"] = default_config
        assembly_config = app.registry.get(session_id)["assembly_config"]
        self.render("config.html",css=css,assembly_config=assembly_config)
    def post(self):
        parts_list,session_id = self.get_parts_list()
        app = webapp2.get_app()
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
                    application.registry[session_id]["assembly_config"][key] = self.request.POST.get('key')
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