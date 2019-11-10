#!/usr/bin/python

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
import logging
import csv

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
#Optionally render a welcome page
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

    try:
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

        #Updates the parts_list by extracting information stored in the HTML DOM, which is modified as the user interacts with the web page
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
    except Exception as e:
        logging.error(e)
        self.redirect("/error")

'''
Development note: the webapp2 documentation specifies that retrieving information from the
web page uses the syntax "self.request.get()". While running in its own server, it seems that
different syntax must be used. Use "self.request.POST.get()" for post method data, and
"self.request.GET.get()" for get method data
'''
try:
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
        primer_forward_tm = 0
        primer_reverse_tm = 0

    #Object for combinatorial assembly; one MultiPart contains multiple parts at equivalent assembly junctions
    class MultiPart:
        def __init__(self,name,parts):
            self.name = name
            self.parts = parts
            self.type = "MultiPart"
        assembly_method = ""

    #Creates a unique session_id for a user to be stored as a cookie
    def make_session_id():
        id_one = {1:'A',2:'B',3:'C',4:'D',5:'E',
        6:'F',7:'G',8:'H',9:'I',10:'J',11:'K',
        12:'L',13:'M',14:'N',15:'O',16:'P',
        17:'Q',18:'R',19:'S',20:'T',
        21:'U',22:'V',23:'W',24:'X',25:'Y',
        26:'Z'
        }
        id_two = {1:'A',2:'B',3:'C',4:'D',5:'E',
        6:'F',7:'G',8:'H',9:'I',10:'J',11:'K',
        12:'L',13:'M',14:'N',15:'O',16:'P',
        17:'Q',18:'R',19:'S',20:'T',
        21:'U',22:'V',23:'W',24:'X',25:'Y',
        26:'Z'
        }
        from random import randrange
        import time
        two = id_two[randrange(1,27)]
        one = id_one[randrange(1,27)]
        number = str(randrange(0,999999999999))
        return one+two+number

    #Takes the current build and assembles it into an XML file in the SBOL data format
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

    #Creates a bridge for LCR in the junction between a sequence and the next sequence in the build
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

    #Creates a builds list from a parts list. Parts lists with MultiPart objects are "unzipped" into multiple builds lists. 
    #For example, a parts list with "A", "B" or "b", and "C" will create a builds list with "A", "B", "C" parts list and "A", "b", "C" parts list
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

    #Optimizes the overhangs used in Golden Gate assembly
    def golden_gate_optimization(parts_list,backbone_sequence,gg_overhangs,bsaI_combs=[]):
        #seq_pairs = {0:"ccct",1:"gctc",2:"cggt",3:"gtgc",4:"agcg",5:"ctgt",6:"tgct",7:"atgg",8:"gact",9:"ggac",10:"tccg",11:"ccag",12:"cagc",13:"gttg",14:"cgaa",15:"ccat"}
        #golden_gate_overhangs = [
        #            "ccct","gctc","cggt","gtgc","agcg","ctgt","tgct","atgg","gact","ggac","tccg","ccag","cagc","gttg","cgaa","ccat"
        #        ]
        golden_gate_overhangs=gg_overhangs
        seq_matches = []
        for x in range(len(parts_list)+1):
            seq_matches.append([])
            if x == 0:
                for overhang in golden_gate_overhangs:
                    if overhang in backbone_sequence[-35:]:
                        seq_matches[x].append(overhang)
                    elif overhang in parts_list[x].sequence[:35]:
                        seq_matches[x].append(overhang)
            elif x == len(parts_list):
                for overhang in golden_gate_overhangs:
                    if overhang in parts_list[x-1].sequence[-35:]:
                        seq_matches[x].append(overhang)
                    elif overhang in backbone_sequence[:35]:
                        seq_matches[x].append(overhang)
            else:
                for overhang in golden_gate_overhangs:
                    if overhang in parts_list[x-1].sequence[-35:]:
                        seq_matches[x].append(overhang)
                    elif overhang in parts_list[x].sequence[:35]:
                        seq_matches[x].append(overhang)
        if bsaI_combs != []:
            for comb in bsaI_combs:
                seq_matches.append(comb)
        combs = []
        for x in itertools.product(*seq_matches):
            combs.append(x)
        for comb in combs:
            if len(comb) == len(set(comb)):
                return comb
        #if there are no possible combinations
        return None

    #Given a sequence and an index where there is a BsaI site (as well as the type of the part), alter the sequence in such a way as to generate a silent mutation to remove the BsaI site
    def silent_mut(sequence,bsa_index,type):
        if type != "CDS":
            if sequence[bsa_index+1] == "g":
                sequence = sequence[:bsa_index] + "gctctc" + sequence[bsa_index+6:]
            else:
                sequence = sequence[:bsa_index] + "gagagc" + sequence[bsa_index+6:]
        else:
            start_index = sequence.find("atg")
            stop_index = 99999999999999
            stop_codon = ""
            for i in range(start_index,len(sequence),3):
                if sequence[i:i+3] == "taa":
                    if i < stop_index:
                        stop_index=i
                        stop_codon="taa"
                elif sequence[i:i+3] == "tag":
                    if i < stop_index:
                        stop_index=i
                        stop_codon="tag"
                elif sequence[i:i+3] == "tga":
                    if i < stop_index:
                        stop_index=i
                        stop_codon="tga"
            if bsa_index < start_index or bsa_index > stop_index:
                if sequence[bsa_index+1] == "g":
                    sequence = sequence[:bsa_index] + "gcactc" + sequence[bsa_index+6:]
                else:
                    sequence = sequence[:bsa_index] + "gagtgc" + sequence[bsa_index+6:]
            elif bsa_index > start_index and bsa_index < stop_index:
                for i in range(start_index,len(sequence),3):
                    if i >= bsa_index:
                        bsa_codon_index=i
                        bsa_codon=sequence[i:i+3]
                        break
                codon_table = [["gct","gcc","gca","gcg"],["cgt","cgc","cga","cgg","aga","agg"],["aat","aac"],["gat","gac"],["tgt","tgc"],
                               ["caa","cag"],["gaa","gag"],["ggt","ggc","gga","ggg"],["cat","cac"],["att","atc","ata"],
                               ["tta","ttg","ctt","ctc","cta","ctg"],["aaa","aag"],["ttt","ttc"],["cct","ccc","cca","ccg"],
                               ["tct","tcc","tca","tcg","agt","agc"],["tat","tac"],["gtt","gtc","gta","gtg"]]
                for aa in codon_table:
                    if bsa_codon in aa:
                        for codon in aa:
                            if codon != bsa_codon:
                                new_codon = codon
                sequence = sequence[:bsa_codon_index] + new_codon + sequence[bsa_codon_index+3:]
        return sequence


    # '/' path
    class MainHandler(Handler):
        def get(self):
            parts_list,session_id = self.get_parts_list()
            self.render("main_page.html",golden_gate_error="",css=css,js=js,parts_list=parts_list)
        def post(self):
            try:
                #If "Delete construct" clicked, deletes existing cookie and assigns a new one
                if self.request.POST.get("delete_map")=="TERMINATE":
                    self.response.delete_cookie("sessionid")
                    self.redirect("/")
    
    			#Extract information from the main page POST request. If some variables on the front or back end have not been initialized, initialize them.
    			#Import file inputs if files are input from the main page
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

                #If "Save Construct" clicked, download the XML document per generateSBOLdoc()
                if self.request.POST.get("save_construct") == "save_construct":
                    parts_list,session_id = self.update_part_list()
                    builds_list = builds(parts_list)
                    generateSBOLdoc(builds_list,session_id)
                    self.redirect("/construct_download")

                #Read in an XML document generated from "Save Construct"
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
    			
    			#Initialize default config settings
                if is_local:
                    for record in SeqIO.parse("templates/pET-26b.fa","fasta"):
                        default_backbone = record
                else:
                    for record in SeqIO.parse("/var/www/ibiocad/iBioCAD/templates/pET-26b.fa","fasta"):
                        default_backbone = record
                default_config = {"backbone":default_backbone,"Golden_gate_method":"scarless_assembly","backbone_primers_tm":[mt.Tm_NN(default_backbone.seq[:40]),mt.Tm_NN(default_backbone[-40:].seq)],"backbone_primers":[default_backbone.seq[:40],default_backbone.seq[-40:]]
                                  ,"Primer_optimization":"both","primer_tm":[52,60]}

                #Read in variables from the config page
                parts_list,session_id = self.get_parts_list()
                application = webapp2.get_app()
                if "assembly_config" not in application.registry.get(session_id).keys() or application.registry[session_id]["assembly_config"] is None:
                    application.registry[session_id]["assembly_config"] = default_config
                assembly_config = application.registry.get(session_id)["assembly_config"]
                backbone_sequence = assembly_config["backbone"].seq
                new_backbone_sequence = backbone_sequence
                backbone_primers_tm = assembly_config["backbone_primers_tm"]
                backbone_primers = assembly_config["backbone_primers"]
                application.registry[session_id]["new_backbone_sequence"] = new_backbone_sequence
                golden_gate_method = assembly_config["Golden_gate_method"]
                primer_optimization = assembly_config["Primer_optimization"]
                primer_tm_range = assembly_config["primer_tm"]

                #If "Assemble" is clicked, run computation based on the chosen Assembly Method

                #Run Yeast Assembly
                #Based on the homology distance, generates primer sequences and modifies the DNA sequences and backbone sequences
                if self.request.POST.get("assembly_method") == "Yeast_Assembly":
                    parts_list,session_id = self.update_part_list()
                    builds_list = builds(parts_list)
                    app = webapp2.get_app()
                    app.registry[session_id]['builds_list'] = builds_list
                    new_backbone_sequence = backbone_sequence
                    new_backbone_primers = [new_backbone_sequence[:40],new_backbone_sequence[-40:]]
                    for unpacked_list in builds_list:
                        backbone_list = []
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
                        new_backbone_sequence = unpacked_list[-1].sequence[-40:] + new_backbone_sequence + unpacked_list[0].sequence[:40]
                        new_backbone_primers = [new_backbone_sequence[:40],new_backbone_sequence[-40:]]
                        backbone_list.append([new_backbone_sequence,new_backbone_primers[0],new_backbone_primers[1]])
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

                    #Optimize primer Tm based on config settings
                    if primer_optimization == "range":
                        for unpacked_list in builds_list:
                            for part in unpacked_list:
                                if mt.Tm_NN(part.primer_forward) < primer_tm_range[0]:
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                else:
                                    for i in range(len(part.primer_forward),15,-1):
                                        if mt.Tm_NN(part.primer_forward[:i]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i]) <= primer_tm_range[1]:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                            part.primer_forward = part.primer_forward[:i]
                                            break
                                        else:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                if mt.Tm_NN(part.primer_reverse) < primer_tm_range[0]:
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                else:
                                    for i in range(len(part.primer_reverse),15,-1):
                                        if mt.Tm_NN(part.primer_reverse[:i]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:i]) <= primer_tm_range[1]:
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:i])
                                            part.primer_reverse = part.primer_reverse[:i]
                                            break
                                        else:
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        breakit=False
                        for x in range(40,15,-1):
                            for z in range(40,15,-1):
                                if mt.Tm_NN(new_backbone_sequence[:x])>=primer_tm_range[0] and mt.Tm_NN(new_backbone_sequence[:x])<=primer_tm_range[1] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]) >= primer_tm_range[0] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]) <= primer_tm_range[1]:
                                    backbone_primers = [new_backbone_sequence[:x],reverse_complement(new_backbone_sequence[-40:])[:z]]
                                    backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:x]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z])]
                                    breakit=True
                                if breakit:
                                    break
                            if breakit:
                                break
                        if not breakit:
                            backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                            backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:40]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:]))]
                    elif primer_optimization == "near":
                        for unpacked_list in builds_list:
                            for part in unpacked_list:
                                breakit=False
                                for i in range(len(part.primer_forward),15,-1):
                                    for j in range(len(part.primer_reverse),15,-1):
                                        if abs(mt.Tm_NN(part.primer_forward[:i])-mt.Tm_NN(part.primer_reverse[:j])) <=5:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                            part.primer_forward = part.primer_forward[:i]
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j])
                                            part.primer_reverse = part.primer_reverse[:j]
                                            breakit=True
                                        if breakit:
                                            break
                                    if breakit:
                                        break
                                if not breakit:
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        breakit=False
                        for x in range(40,15,-1):
                            for z in range(40,15,-1):
                                if abs(mt.Tm_NN(new_backbone_sequence[:x])-mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]))<=5:
                                    backbone_primers = [new_backbone_sequence[:x],reverse_complement(new_backbone_sequence[-40:])[:z]]
                                    backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:x]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z])]
                                    breakit=True
                                if breakit:
                                    break
                            if breakit:
                                break
                        if not breakit:
                            backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                            backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:40]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:]))]
                    elif primer_optimization == "both":
                        for unpacked_list in builds_list:
                            for part in unpacked_list:
                                if mt.Tm_NN(part.primer_forward) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse) < primer_tm_range[0]:
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                elif mt.Tm_NN(part.primer_forward) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse) >= primer_tm_range[0]:
                                    breakit=False
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                    for j in range(len(part.primer_reverse),15,-1):
                                        if abs(part.primer_forward_tm-mt.Tm_NN(part.primer_reverse[:j])) <=5:
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j])
                                            part.primer_reverse = part.primer_reverse[:j]
                                            breakit=True
                                    if not breakit:
                                        part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                elif mt.Tm_NN(part.primer_reverse) < primer_tm_range[0] and mt.Tm_NN(part.primer_forward) >= primer_tm_range[0]:
                                    breakit=False
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                    for i in range(len(part.primer_forward),15,-1):
                                        if abs(part.primer_reverse_tm-mt.Tm_NN(part.primer_forward[:i])) <=5:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                            part.primer_forward = part.primer_forward[:i]
                                            breakit=True
                                    if not breakit:
                                        part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                else:
                                    breakit=False
                                    for i in range(len(part.primer_forward),15,-1):
                                        for j in range(len(part.primer_reverse),15,-1):
                                            if abs(mt.Tm_NN(part.primer_forward[:i])-mt.Tm_NN(part.primer_reverse[:j]))<=5 and mt.Tm_NN(part.primer_forward[:i]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i]) <= primer_tm_range[1] and mt.Tm_NN(part.primer_reverse[:j]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:j]) <= primer_tm_range[1]:
                                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                                part.primer_forward = part.primer_forward[:i]
                                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j])
                                                part.primer_reverse = part.primer_reverse[:j]
                                                breakit=True
                                            if breakit:
                                                break
                                        if breakit:
                                            break
                                    if not breakit:
                                        part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                        part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        breakit=False
                        for x in range(40,15,-1):
                            for z in range(40,15,-1):
                                if abs(mt.Tm_NN(new_backbone_sequence[:x])-mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]))<=5 and mt.Tm_NN(new_backbone_sequence[:x])>=primer_tm_range[0] and mt.Tm_NN(new_backbone_sequence[:x])<=primer_tm_range[1] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]) >= primer_tm_range[0] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]) <= primer_tm_range[1]:
                                    backbone_primers = [new_backbone_sequence[:x],reverse_complement(new_backbone_sequence[-40:])[:z]]
                                    backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:x]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z])]
                                    breakit=True
                                if breakit:
                                    break
                            if breakit:
                                break
                        if not breakit:
                            backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                            backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:40]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:]))]
                    else:
                        for unpacked_list in builds_list:
                            for part in unpacked_list:
                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                        backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:40]),reverse_complement(new_backbone_sequence[-40:])]

                    #Register new parts and direct to output page
                    app.registry[session_id]["backbone_primers"] = backbone_primers
                    app.registry[session_id]["backbone_primers_tm"] = backbone_primers_tm
                    parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
                    app.registry[session_id]['builds_list'] = builds_list
                    app.registry[session_id]["backbone_list"] = backbone_list
                    app.registry[session_id]["new_backbone_sequence"] = new_backbone_sequence
                    self.redirect("/assembly")

                #Run Gibson Assembly
                #Based on the homology distance, generates primer sequences and modifies the DNA sequences and backbone sequences
                if self.request.POST.get("assembly_method") == "Gibson_Assembly":
                    parts_list,session_id = self.update_part_list()
                    builds_list = builds(parts_list)
                    app = webapp2.get_app()
                    app.registry[session_id]['builds_list'] = builds_list
                    new_backbone_sequence = backbone_sequence
                    new_backbone_primers = [new_backbone_sequence[:40],new_backbone_sequence[-40:]]
                    for unpacked_list in builds_list:
                        backbone_list = []
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
                        new_backbone_sequence = unpacked_list[-1].sequence[-25:] + new_backbone_sequence + unpacked_list[0].sequence[:25]
                        new_backbone_primers = [new_backbone_sequence[:40],new_backbone_sequence[-40:]]
                        backbone_list.append([new_backbone_sequence,new_backbone_primers[0],new_backbone_primers[1]])
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

                    #Optimize primer Tm based off the config settings
                    if primer_optimization == "range":
                        for unpacked_list in builds_list:
                            for part in unpacked_list:
                                if mt.Tm_NN(part.primer_forward) < primer_tm_range[0]:
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                else:
                                    for i in range(len(part.primer_forward),15,-1):
                                        if mt.Tm_NN(part.primer_forward[:i]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i]) <= primer_tm_range[1]:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                            part.primer_forward = part.primer_forward[:i]
                                            break
                                        else:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                if mt.Tm_NN(part.primer_reverse) < primer_tm_range[0]:
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                else:
                                    for i in range(len(part.primer_reverse),15,-1):
                                        if mt.Tm_NN(part.primer_reverse[:i]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:i]) <= primer_tm_range[1]:
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:i])
                                            part.primer_reverse = part.primer_reverse[:i]
                                            break
                                        else:
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        breakit=False
                        for x in range(40,15,-1):
                            for z in range(40,15,-1):
                                if mt.Tm_NN(new_backbone_sequence[:x])>=primer_tm_range[0] and mt.Tm_NN(new_backbone_sequence[:x])<=primer_tm_range[1] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]) >= primer_tm_range[0] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]) <= primer_tm_range[1]:
                                    backbone_primers = [new_backbone_sequence[:x],reverse_complement(new_backbone_sequence[-40:])[:z]]
                                    backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:x]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z])]
                                    breakit=True
                                if breakit:
                                    break
                            if breakit:
                                break
                        if not breakit:
                            backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                            backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:40]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:]))]
                    elif primer_optimization == "near":
                        for unpacked_list in builds_list:
                            for part in unpacked_list:
                                breakit=False
                                for i in range(len(part.primer_forward),15,-1):
                                    for j in range(len(part.primer_reverse),15,-1):
                                        if abs(mt.Tm_NN(part.primer_forward[:i])-mt.Tm_NN(part.primer_reverse[:j])) <=5:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                            part.primer_forward = part.primer_forward[:i]
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j])
                                            part.primer_reverse = part.primer_reverse[:j]
                                            breakit=True
                                        if breakit:
                                            break
                                    if breakit:
                                        break
                                if not breakit:
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        breakit=False
                        for x in range(40,15,-1):
                            for z in range(40,15,-1):
                                if abs(mt.Tm_NN(new_backbone_sequence[:x])-mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]))<=5:
                                    backbone_primers = [new_backbone_sequence[:x],reverse_complement(new_backbone_sequence[-40:])[:z]]
                                    backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:x]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z])]
                                    breakit=True
                                if breakit:
                                    break
                            if breakit:
                                break
                        if not breakit:
                            backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                            backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:40]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:]))]
                    elif primer_optimization == "both":
                        for unpacked_list in builds_list:
                            for part in unpacked_list:
                                if mt.Tm_NN(part.primer_forward) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse) < primer_tm_range[0]:
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                elif mt.Tm_NN(part.primer_forward) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse) >= primer_tm_range[0]:
                                    breakit=False
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                    for j in range(len(part.primer_reverse),15,-1):
                                        if abs(part.primer_forward_tm-mt.Tm_NN(part.primer_reverse[:j])) <=5:
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j])
                                            part.primer_reverse = part.primer_reverse[:j]
                                            breakit=True
                                    if not breakit:
                                        part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                elif mt.Tm_NN(part.primer_reverse) < primer_tm_range[0] and mt.Tm_NN(part.primer_forward) >= primer_tm_range[0]:
                                    breakit=False
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                    for i in range(len(part.primer_forward),15,-1):
                                        if abs(part.primer_reverse_tm-mt.Tm_NN(part.primer_forward[:i])) <=5:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                            part.primer_forward = part.primer_forward[:i]
                                            breakit=True
                                    if not breakit:
                                        part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                else:
                                    breakit=False
                                    for i in range(len(part.primer_forward),15,-1):
                                        for j in range(len(part.primer_reverse),15,-1):
                                            if abs(mt.Tm_NN(part.primer_forward[:i])-mt.Tm_NN(part.primer_reverse[:j]))<=5 and mt.Tm_NN(part.primer_forward[:i]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i]) <= primer_tm_range[1] and mt.Tm_NN(part.primer_reverse[:j]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:j]) <= primer_tm_range[1]:
                                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                                part.primer_forward = part.primer_forward[:i]
                                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j])
                                                part.primer_reverse = part.primer_reverse[:j]
                                                breakit=True
                                            if breakit:
                                                break
                                        if breakit:
                                            break
                                    if not breakit:
                                        part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                        part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        breakit=False
                        for x in range(40,15,-1):
                            for z in range(40,15,-1):
                                if abs(mt.Tm_NN(new_backbone_sequence[:x])-mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]))<=5 and mt.Tm_NN(new_backbone_sequence[:x])>=primer_tm_range[0] and mt.Tm_NN(new_backbone_sequence[:x])<=primer_tm_range[1] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]) >= primer_tm_range[0] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z]) <= primer_tm_range[1]:
                                    backbone_primers = [new_backbone_sequence[:x],reverse_complement(new_backbone_sequence[-40:])[:z]]
                                    backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:x]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:])[:z])]
                                    breakit=True
                                if breakit:
                                    break
                            if breakit:
                                break
                        if not breakit:
                            backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                            backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:40]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-40:]))]
                    else:
                        for unpacked_list in builds_list:
                            for part in unpacked_list:
                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                        backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:40]),reverse_complement(new_backbone_sequence[-40:])]

                    #Register new parts and direct to output page
                    app.registry[session_id]["backbone_primers"] = backbone_primers
                    app.registry[session_id]["backbone_primers_tm"] = backbone_primers_tm
                    parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
                    app.registry[session_id]['builds_list'] = builds_list
                    app.registry[session_id]["backbone_list"] = backbone_list
                    app.registry[session_id]["new_backbone_sequence"] = new_backbone_sequence
                    self.redirect("/assembly")

                #Run Ligase Cycling Reaction
                #Designs LCR bridges based off the LCR method (de Kok et al, 2014)
                if self.request.POST.get("assembly_method") == "LCR":
                    parts_list,session_id = self.update_part_list()
                    builds_list = builds(parts_list)
                    application = webapp2.get_app()
                    application.registry[session_id]['builds_list'] = builds_list
                    backbone_list=[]
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
                    
                    #Register new parts and direct to output page
                    application.registry[session_id]["backbone_primers"] = backbone_primers
                    application.registry[session_id]["backbone_primers_tm"] = backbone_primers_tm
                    parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
                    application.registry[session_id]['builds_list'] = builds_list
                    application.registry[session_id]["backbone_list"] = backbone_list
                    self.redirect("/assembly")

                #Run Golden Gate Assembly Design
                #Based off the settings, particularly regular or scarless Golden Gate assembly and whether to remove BsaI sequences internal to the parts, the computation will be different
                golden_gate_error=""
                if self.request.POST.get("assembly_method") == "Type_II_Restriction_Enzyme":
                    remove_bsaI = self.request.POST.get("remove_bsaI")
                    parts_list,session_id = self.update_part_list()
                    builds_list = builds(parts_list)
                    app = webapp2.get_app()
                    app.registry[session_id]['builds_list'] = builds_list
                    new_backbone_sequence = backbone_sequence
                    new_backbone_primers = [new_backbone_sequence[:40],new_backbone_sequence[-40:]]

                    #Run regular Golden Gate assembly design
                    if golden_gate_method == "regular_assembly":
                        new_backbone_sequence = app.registry.get(session_id)["new_backbone_sequence"]
                        new_builds_list = []
                        for unpacked_list in builds_list:
                            new_unpacked_list = []
                            contains_bsaI = False
                            #Check for BsaI sites in the sequences
                            for part in unpacked_list:
                                part.assembly_method = "Type_II_Restriction_Enzyme"
                                part.sequence = part.sequence.lower()
                                if "ggtctc" in part.sequence or "gagacc" in part.sequence or "ggtctc" in backbone_sequence or "gagacc" in backbone_sequence:
                                    golden_gate_error = "BsaI_in_seq" +"|"+ part.name
                                    contains_bsaI = True
                            #Perform computation assuming BsaI sites must be removed
                            if remove_bsaI == "yes":
                            	#Reserve overhangs for removing internal BsaI sites based on a set of GG linkers
                                def gg_reg_opt(unpacked_list,gg_overhangs):
                                    new_backbone_sequence = app.registry.get(session_id)["new_backbone_sequence"]
                                    golden_gate_overhangs = gg_overhangs
                                    golden_gate_error = ""
                                    reserved_overhangs = []
                                    bsaI_removal_map = {}
                                    bsaI_removal_combs = []
                                    bsaI_removal_key = []
                                    for part in unpacked_list:
                                        part.sequence = part.sequence.lower()
                                        bsa_index = []
                                        ii=0
                                        while len(bsa_index) != (part.sequence.count("ggtctc")+part.sequence.count("gagacc")) and ii<10*len(unpacked_list):
                                            if "ggtctc" in part.sequence and "gagacc" not in part.sequence:
                                                if bsa_index == []:
                                                    bsa_index.append(part.sequence.find("ggtctc"))
                                                else:
                                                    bsa_index.append(part.sequence.find("ggtctc",bsa_index[-1]+1))
                                            elif "gagacc" in part.sequence and "ggtctc" not in part.sequence:
                                                if bsa_index == []:
                                                    bsa_index.append(part.sequence.find("gagacc"))
                                                else:
                                                    bsa_index.append(part.sequence.find("gagacc",bsa_index[-1]+1))
                                            elif "ggtctc" in part.sequence and "gagacc" in part.sequence:
                                                if bsa_index == []:
                                                    bsa_index.append(min(part.sequence.find("ggtctc"),part.sequence.find("gagacc")))
                                                else:
                                                    bsa_index.append(min(part.sequence.find("ggtctc",bsa_index[-1]+1),part.sequence.find("gagacc",bsa_index[-1]+1)))
                                            ii+=1
                                        if "ggtctc" in part.sequence or "gagacc" in part.sequence:
                                            for xx in bsa_index:
                                                bsaI_removal_key.append(part)
                                            for i in range(len(bsa_index)):
                                                bsaI_removal_combs.append([])
                                                for overhang in golden_gate_overhangs:
                                                    if overhang in part.sequence[bsa_index[i]-30:bsa_index[i]+37] and ((part.sequence[bsa_index[i]-30:bsa_index[i]+37].find(overhang)+bsa_index[i]-30)<(bsa_index[i]-3) or (part.sequence[bsa_index[i]-30:bsa_index[i]+37].find(overhang)+bsa_index[i]-30)>=(bsa_index[i]+6)):
                                                        bsaI_removal_combs[-1].append(overhang)
                                    new_backbone_sequence = new_backbone_sequence.lower()
                                    bsa_index = []
                                    ii = 0
                                    while len(bsa_index) != (backbone_sequence.count("ggtctc")+backbone_sequence.count("gagacc")) and ii<10:
                                        if "ggtctc" in backbone_sequence and "gagacc" not in backbone_sequence:
                                            if bsa_index == []:
                                                bsa_index.append(backbone_sequence.find("ggtctc"))
                                            else:
                                                bsa_index.append(backbone_sequence.find("ggtctc",bsa_index[-1]+1))
                                        elif "gagacc" in backbone_sequence and "ggtctc" not in backbone_sequence:
                                            if bsa_index == []:
                                                bsa_index.append(backbone_sequence.find("gagacc"))
                                            else:
                                                bsa_index.append(backbone_sequence.find("gagacc",bsa_index[-1]+1))
                                        elif "ggtctc" in backbone_sequence and "gagacc" in backbone_sequence:
                                            if bsa_index == []:
                                                bsa_index.append(min(backbone_sequence.find("ggtctc"),backbone_sequence.find("gagacc")))
                                            else:
                                                if backbone_sequence.find("ggtctc",bsa_index[-1]+1) == -1 or backbone_sequence.find("gagacc",bsa_index[-1]+1) == -1:
                                                    bsa_index.append(max(backbone_sequence.find("ggtctc",bsa_index[-1]+1),backbone_sequence.find("gagacc",bsa_index[-1]+1)))
                                                else:
                                                    bsa_index.append(min(backbone_sequence.find("ggtctc",bsa_index[-1]+1),backbone_sequence.find("gagacc",bsa_index[-1]+1)))
                                        ii+=1
                                    if "ggtctc" in backbone_sequence or "gagacc" in backbone_sequence:
                                        for xx in bsa_index:
                                            bsaI_removal_key.append("backbone")
                                        for i in range(len(bsa_index)):
                                            bsaI_removal_combs.append([])
                                            for overhang in golden_gate_overhangs:
                                                if overhang in backbone_sequence[bsa_index[i]-30:bsa_index[i]+37] and ((backbone_sequence[bsa_index[i]-30:bsa_index[i]+37].find(overhang)+bsa_index[i]-30)<(bsa_index[i]-3) or (backbone_sequence[bsa_index[i]-30:bsa_index[i]+37].find(overhang)+bsa_index[i]-30)>=(bsa_index[i]+6)):
                                                    bsaI_removal_combs[-1].append(overhang)
                                    combs = []
                                    for x in itertools.product(*bsaI_removal_combs):
                                        combs.append(x)
                                    for comb in combs:
                                        if len(comb) == len(set(comb)):
                                            reserved_overhangs = comb
                                    if reserved_overhangs != []:
                                        for i in range(len(bsaI_removal_key)):
                                            if bsaI_removal_key[i]=="backbone" and "backbone" in bsaI_removal_map.keys():
                                                bsaI_removal_map["backbone"].append(reserved_overhangs[i])
                                            elif bsaI_removal_key[i]=="backbone" and "backbone" not in bsaI_removal_map.keys():
                                                bsaI_removal_map["backbone"] = [reserved_overhangs[i]]
                                            elif bsaI_removal_key[i].name in bsaI_removal_map.keys():
                                               bsaI_removal_map[bsaI_removal_key[i].name].append(reserved_overhangs[i])
                                            else:
                                                bsaI_removal_map[bsaI_removal_key[i].name] = [reserved_overhangs[i]]
                                    else:
                                        golden_gate_error = "no_efficient_overhang_combinations"
                                        bsaI_removal_key = []
                                    return reserved_overhangs,bsaI_removal_map,golden_gate_error
                                #Try multiple sets of overhangs until a suitable set has been found
                                breakit = False
                                for i in range(10,51):
                                    for j in range(5):
                                        with open('/var/www/ibiocad/iBioCAD/overhangsets/setsof%s.csv'%i,'r') as f:
                                            reader = csv.reader(f, delimiter=",")
                                            temp = list(reader)[1:]
                                        gg_overhangs = []
                                        for x in range(i):
                                            gg_overhangs.append(temp[j][x].lower())
                                        if gg_reg_opt(unpacked_list,gg_overhangs)[2] == "" and (len(gg_reg_opt(unpacked_list,gg_overhangs)[0])+len(unpacked_list))<i :
                                            breakit = True
                                            reserved_overhangs,bsaI_removal_map,golden_gate_error = gg_reg_opt(unpacked_list,gg_overhangs)
                                        if breakit:
                                            break
                                    if breakit:
                                        break
                            else:
                                bsaI_removal_map = {}
                                oh_num = max(len(unpacked_list),10)
                                with open('/var/www/ibiocad/iBioCAD/overhangsets/setsof%s.csv'%oh_num,'r') as f:
                                    reader = csv.reader(f, delimiter=",")
                                    temp = list(reader)[1:]
                                gg_overhangs = []
                                for x in range(oh_num):
                                    gg_overhangs.append(temp[0][x].lower())
                            golden_gate_overhangs = gg_overhangs
                            #Design flanking elements for Golden Gate assembly and primer design
                            for i in range(len(unpacked_list)):
                                if golden_gate_error != "" or remove_bsaI == "no":
                                    break
                                if remove_bsaI == "yes":
                                    usable_overhangs = list(set(golden_gate_overhangs)-set(reserved_overhangs))
                                else:
                                    usable_overhangs = list(set(golden_gate_overhangs))
                                if len(unpacked_list) > len(usable_overhangs):
                                    break
                                unpacked_list[i].primer_forward = "aaggtctca" + usable_overhangs[i] + unpacked_list[i].sequence[:20]
                                unpacked_list[i].primer_reverse = reverse_complement(unpacked_list[i].sequence[-20:] + usable_overhangs[i+1] + "agagaccaa")
                                unpacked_list[i].sequence = "aaggtctca" + usable_overhangs[i] + unpacked_list[i].sequence + usable_overhangs[i+1] + "agagaccaa"
                                #Remove internal BsaI sites from non-backbone parts
                                if unpacked_list[i].name in bsaI_removal_map.keys():
                                    bsa_index = []
                                    iii=0
                                    while len(bsa_index) != (unpacked_list[i].sequence.count("ggtctc")+unpacked_list[i].sequence.count("gagacc"))-2 and iii<10*len(unpacked_list):
                                        if "ggtctc" in unpacked_list[i].sequence[13:] and "gagacc" not in unpacked_list[i].sequence[:len(unpacked_list[i].sequence)-13]:
                                            if bsa_index == []:
                                                bsa_index.append(unpacked_list[i].sequence.find("ggtctc",13))
                                            else:
                                                bsa_index.append(unpacked_list[i].sequence.find("ggtctc",bsa_index[-1]+1))
                                        elif "gagacc" in unpacked_list[i].sequence[:len(unpacked_list[i].sequence)-13] and "ggtctc" not in unpacked_list[i].sequence[13:]:
                                            if bsa_index == []:
                                                bsa_index.append(unpacked_list[i].sequence.find("gagacc",0,len(unpacked_list[i].sequence)-13))
                                            else:
                                                bsa_index.append(unpacked_list[i].sequence.find("gagacc",bsa_index[-1]+1))
                                        elif "ggtctc" in unpacked_list[i].sequence[13:] and "gagacc" in unpacked_list[i].sequence[:len(unpacked_list[i].sequence)-13]:
                                            if bsa_index == []:
                                                bsa_index.append(min(unpacked_list[i].sequence.find("ggtctc",13),unpacked_list[i].sequence.find("gagacc",0,len(unpacked_list[i].sequence)-13)))
                                            else:
                                                bsa_index.append(min(unpacked_list[i].sequence.find("ggtctc",bsa_index[-1]+1),unpacked_list[i].sequence.find("gagacc",bsa_index[-1]+1)))
                                        iii+=1
                                    for x in range(len(bsaI_removal_map[unpacked_list[i].name])):
                                        overhang = bsaI_removal_map[unpacked_list[i].name][x]
                                        if x == 0:
                                            if bsa_index[x]-30 < 0:
                                                if overhang in unpacked_list[i].sequence[0:bsa_index[x]]:
                                                    overhang_index = unpacked_list[i].sequence.find(overhang,0,bsa_index[x])
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"a",unpacked_list[i].type,unpacked_list[i].sequence[:overhang_index]+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = unpacked_list[i].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[i].sequence[overhang_index-20:overhang_index]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"b",unpacked_list[i].type,"aaggtctca"+silent_mut(unpacked_list[i].sequence[overhang_index:],bsa_index[x]-overhang_index,unpacked_list[i].type)))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + silent_mut(unpacked_list[i].sequence[overhang_index:],bsa_index[x]-overhang_index,unpacked_list[i].type)[:24]
                                                    new_unpacked_list[-1].primer_reverse = unpacked_list[i].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                elif overhang in unpacked_list[i].sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                    overhang_index = unpacked_list[i].sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"a",unpacked_list[i].type,silent_mut(unpacked_list[i].sequence[:overhang_index],bsa_index[x],unpacked_list[i].type)+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = unpacked_list[i].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(silent_mut(unpacked_list[i].sequence[:overhang_index],bsa_index[x],unpacked_list[i].type)[-20:]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"b",unpacked_list[i].type,"aaggtctca"+unpacked_list[i].sequence[overhang_index:]))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[i].sequence[overhang_index:overhang_index+24]
                                                    new_unpacked_list[-1].primer_reverse = unpacked_list[i].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                else:
                                                    pass
                                            else:
                                                if overhang in unpacked_list[i].sequence[bsa_index[x]-30:bsa_index[x]]:
                                                    overhang_index = unpacked_list[i].sequence.find(overhang,bsa_index[x]-30,bsa_index[x])
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"a",unpacked_list[i].type,unpacked_list[i].sequence[:overhang_index]+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = unpacked_list[i].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[i].sequence[overhang_index-20:overhang_index]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"b",unpacked_list[i].type,"aaggtctca"+silent_mut(unpacked_list[i].sequence[overhang_index:],bsa_index[x]-overhang_index,unpacked_list[i].type)))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + silent_mut(unpacked_list[i].sequence[overhang_index:],bsa_index[x]-overhang_index,unpacked_list[i].type)[:24]
                                                    new_unpacked_list[-1].primer_reverse = unpacked_list[i].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                elif overhang in unpacked_list[i].sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                    overhang_index = unpacked_list[i].sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"a",unpacked_list[i].type,silent_mut(unpacked_list[i].sequence[:overhang_index],bsa_index[x],unpacked_list[i].type)+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = unpacked_list[i].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(silent_mut(unpacked_list[i].sequence[:overhang_index],bsa_index[x],unpacked_list[i].type)[-20:]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"b",unpacked_list[i].type,"aaggtctca"+unpacked_list[i].sequence[overhang_index:]))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[i].sequence[overhang_index:overhang_index+24]
                                                    new_unpacked_list[-1].primer_reverse = unpacked_list[i].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                else:
                                                    pass
                                        else:
                                            if bsa_index[x]-30 < 0:
                                                if overhang in new_unpacked_list[-1].sequence[0:bsa_index[x]]:
                                                    overhang_index = new_unpacked_list[-1].sequence.find(overhang,0,bsa_index[x])
                                                    new_unpacked_list.append(Part(new_unpacked_list[-1].name+"-"+str(x+1)+"a",new_unpacked_list[-1].type,new_unpacked_list[-1].sequence[:overhang_index]+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = new_unpacked_list[-2].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(new_unpacked_list[-2].sequence[overhang_index-20:overhang_index]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(new_unpacked_list[-2].name+"-"+str(x+1)+"b",new_unpacked_list[-2].type,"aaggtctca"+silent_mut(new_unpacked_list[-2].sequence[overhang_index:],bsa_index[x]-overhang_index,new_unpacked_list[-2].type)))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + silent_mut(new_unpacked_list[-3].sequence[overhang_index:],bsa_index[x]-overhang_index,new_unpacked_list[-3].type)[:24]
                                                    new_unpacked_list[-1].primer_reverse = new_unpacked_list[-3].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    del new_unpacked_list[-3]
                                                elif overhang in new_unpacked_list[-1].sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                    overhang_index = new_unpacked_list[-1].sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                    new_unpacked_list.append(Part(new_unpacked_list[-1].name+"-"+str(x+1)+"a",new_unpacked_list[-1].type,silent_mut(new_unpacked_list[-1].sequence[:overhang_index],bsa_index[x],new_unpacked_list[-1].type)+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = new_unpacked_list[-2].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(silent_mut(new_unpacked_list[-2].sequence[:overhang_index],bsa_index[x],new_unpacked_list[-2].type)[-20:]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(new_unpacked_list[-2].name+"-"+str(x+1)+"b",new_unpacked_list[-2].type,"aaggtctca"+new_unpacked_list[-2].sequence[overhang_index:]))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + new_unpacked_list[-3].sequence[overhang_index:overhang_index+24]
                                                    new_unpacked_list[-1].primer_reverse = new_unpacked_list[-3].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    del new_unpacked_list[-3]
                                                else:
                                                    pass
                                            else:
                                                if overhang in new_unpacked_list[-1].sequence[bsa_index[x]-30:bsa_index[x]]:
                                                    overhang_index = new_unpacked_list[-1].sequence.find(overhang,bsa_index[x]-30,bsa_index[x])
                                                    new_unpacked_list.append(Part(new_unpacked_list[-1].name+"-"+str(x+1)+"a",new_unpacked_list[-1].type,new_unpacked_list[-1].sequence[:overhang_index]+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = new_unpacked_list[-2].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(new_unpacked_list[-2].sequence[overhang_index-20:overhang_index]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(new_unpacked_list[-2].name+"-"+str(x+1)+"b",new_unpacked_list[-2].type,"aaggtctca"+silent_mut(new_unpacked_list[-2].sequence[overhang_index:],bsa_index[x]-overhang_index,new_unpacked_list[-2].type)))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + silent_mut(new_unpacked_list[-3].sequence[overhang_index:],bsa_index[x]-overhang_index,new_unpacked_list[-3].type)[:24]
                                                    new_unpacked_list[-1].primer_reverse = new_unpacked_list[-3].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    del new_unpacked_list[-3]
                                                elif overhang in new_unpacked_list[-1].sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                    overhang_index = new_unpacked_list[-1].sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                    new_unpacked_list.append(Part(new_unpacked_list[-1].name+"-"+str(x+1)+"a",new_unpacked_list[-1].type,silent_mut(new_unpacked_list[-1].sequence[:overhang_index],bsa_index[x],new_unpacked_list[-1].type)+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = new_unpacked_list[-2].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(silent_mut(new_unpacked_list[-2].sequence[:overhang_index],bsa_index[x],new_unpacked_list[-2].type)[-20:]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(new_unpacked_list[-2].name+"-"+str(x+1)+"b",new_unpacked_list[-2].type,"aaggtctca"+new_unpacked_list[-2].sequence[overhang_index:]))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + new_unpacked_list[-3].sequence[overhang_index:overhang_index+24]
                                                    new_unpacked_list[-1].primer_reverse = new_unpacked_list[-3].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    del new_unpacked_list[-3]
                                                else:
                                                    pass
                                else:
                                    new_unpacked_list.append(unpacked_list[i])
                            #Remove BsaI sites from backbone sequence
                            if "backbone" in bsaI_removal_map.keys():
                                bsa_index = []
                                ii=0
                                while len(bsa_index) != (new_backbone_sequence.count("ggtctc")+new_backbone_sequence.count("gagacc"))-2 and ii<10:
                                    if "ggtctc" in new_backbone_sequence[13:] and "gagacc" not in new_backbone_sequence[:len(new_backbone_sequence)-13]:
                                        if bsa_index == []:
                                            bsa_index.append(new_backbone_sequence.find("ggtctc",13))
                                        else:
                                            bsa_index.append(new_backbone_sequence.find("ggtctc",bsa_index[-1]+1))
                                    elif "gagacc" in new_backbone_sequence[:len(new_backbone_sequence)-13] and "ggtctc" not in new_backbone_sequence[13:]:
                                        if bsa_index == []:
                                            bsa_index.append(new_backbone_sequence.find("gagacc",0,len(new_backbone_sequence)-13))
                                        else:
                                            bsa_index.append(new_backbone_sequence.find("gagacc",bsa_index[-1]+1))
                                    elif "ggtctc" in new_backbone_sequence[13:] and "gagacc" in new_backbone_sequence[:len(new_backbone_sequence)-13]:
                                        if bsa_index == []:
                                            bsa_index.append(min(new_backbone_sequence.find("ggtctc",13),new_backbone_sequence.find("gagacc",0,len(new_backbone_sequence)-13)))
                                        else:
                                            bsa_index.append(min(new_backbone_sequence.find("ggtctc",bsa_index[-1]+1),new_backbone_sequence.find("gagacc",bsa_index[-1]+1)))
                                    ii+=1
                                backbone_list = []
                                new_backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                                for x in range(len(bsaI_removal_map["backbone"])):
                                    backbone_list.append([])
                                    backbone_list.append([])
                                    overhang = bsaI_removal_map("backbone")[x]
                                    if x == 0:
                                        if bsa_index[x]-30 < 0:
                                            if overhang in new_backbone_sequence[0:bsa_index[x]]:
                                                overhang_index = new_backbone_sequence.find(overhang,0,bsa_index[x])
                                                backbone_list[-2].append(new_backbone_sequence[:overhang_index]+overhang+"agagaccaa")
                                                backbone_list[-2].append(new_backbone_primers[0])
                                                backbone_list[-2].append(reverse_complement(new_backbone_sequence[overhang_index-40:overhang_index]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+silent_mut(new_backbone_sequence[overhang_index:],bsa_index[x]-overhang_index,"CDS"))
                                                backbone_list[-1].append("aaggtctca" + silent_mut(new_backbone_sequence[overhang_index:],bsa_index[x]-overhang_index,"CDS")[:44])
                                                backbone_list[-1].append(new_backbone_primers[1])
                                            elif overhang in new_backbone_sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                overhang_index = new_backbone_sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                backbone_list[-2].append(silent_mut(new_backbone_sequence[:overhang_index],bsa_index[x],"CDS")+overhang+"agagaccaa")
                                                backbone_list[-2].append(new_backbone_primers[0])
                                                backbone_list[-2].append(reverse_complement(silent_mut(new_backbone_sequence[:overhang_index],bsa_index[x],"CDS")[-40:]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+new_backbone_sequence[overhang_index:])
                                                backbone_list[-1].append("aaggtctca" + new_backbone_sequence[overhang_index:overhang_index+44])
                                                backbone_list[-1].append(new_backbone_primers[1])
                                            else:
                                                pass
                                        else:
                                            if overhang in new_backbone_sequence[bsa_index[x]-30:bsa_index[x]]:
                                                overhang_index = new_backbone_sequence.find(overhang,bsa_index[x]-30,bsa_index[x])
                                                backbone_list[-2].append(new_backbone_sequence[:overhang_index]+overhang+"agagaccaa")
                                                backbone_list[-2].append(new_backbone_primers[0])
                                                backbone_list[-2].append(reverse_complement(new_backbone_sequence[overhang_index-40:overhang_index]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+silent_mut(new_backbone_sequence[overhang_index:],bsa_index[x]-overhang_index,"CDS"))
                                                backbone_list[-1].append("aaggtctca" + silent_mut(new_backbone_sequence[overhang_index:],bsa_index[x]-overhang_index,"CDS")[:44])
                                                backbone_list[-1].append(new_backbone_primers[1])
                                            elif overhang in new_backbone_sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                overhang_index = new_backbone_sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                backbone_list[-2].append(silent_mut(new_backbone_sequence[:overhang_index],bsa_index[x],"CDS")+overhang+"agagaccaa")
                                                backbone_list[-2].append(new_backbone_primers[0])
                                                backbone_list[-2].append(reverse_complement(silent_mut(new_backbone_sequence[:overhang_index],bsa_index[x],"CDS")[-40:]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+new_backbone_sequence[overhang_index:])
                                                backbone_list[-1].append("aaggtctca" + new_backbone_sequence[overhang_index:overhang_index+44])
                                                backbone_list[-1].append(new_backbone_primers[1])
                                            else:
                                                pass
                                    else:
                                        if bsa_index[x]-30 < 0:
                                            if overhang in backbone_list[-3][0][0:bsa_index[x]]:
                                                overhang_index = backbone_list[-3][0].find(overhang,0,bsa_index[x])
                                                backbone_list[-2].append(backbone_list[-3][0][:overhang_index]+overhang+"agagaccaa")
                                                backbone_list[-2].append(backbone_list[-3][1])
                                                backbone_list[-2].append(reverse_complement(backbone_list[-3][0][overhang_index-40:overhang_index]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+silent_mut(backbone_list[-3][0][overhang_index:],bsa_index[x]-overhang_index,"CDS"))
                                                backbone_list[-1].append("aaggtctca" + silent_mut(backbone_list[-3][0][overhang_index:],bsa_index[x]-overhang_index,"CDS")[:44])
                                                backbone_list[-1].append(backbone_list[-3][2])
                                                del backbone_list[-3]
                                            elif overhang in backbone_list[-3][0][bsa_index[x]+6:bsa_index[x]+37]:
                                                overhang_index = backbone_list[-3][0].find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                backbone_list[-2].append(silent_mut(backbone_list[-3][0][:overhang_index],bsa_index[x],"CDS")+overhang+"agagaccaa")
                                                backbone_list[-2].append(backbone_list[-3][1])
                                                backbone_list[-2].append(reverse_complement(silent_mut(backbone_list[-3][0][:overhang_index],bsa_index[x],"CDS")[-40:]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+backbone_list[-3][0][overhang_index:])
                                                backbone_list[-1].append("aaggtctca" + backbone_list[-3][0][overhang_index:overhang_index+44])
                                                backbone_list[-1].append(backbone_list[-3][2])
                                                del backbone_list[-3]
                                            else:
                                                pass
                                        else:
                                            if overhang in backbone_list[-3][0][bsa_index[x]-30:bsa_index[x]]:
                                                overhang_index = backbone_list[-3][0].find(overhang,bsa_index[x]-30,bsa_index[x])
                                                backbone_list[-2].append(backbone_list[-3][0][:overhang_index]+overhang+"agagaccaa")
                                                backbone_list[-2].append(backbone_list[-3][1])
                                                backbone_list[-2].append(reverse_complement(backbone_list[-3][0][overhang_index-40:overhang_index]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+silent_mut(backbone_list[-3][0][overhang_index:],bsa_index[x]-overhang_index,"CDS"))
                                                backbone_list[-1].append("aaggtctca" + silent_mut(backbone_list[-3][0][overhang_index:],bsa_index[x]-overhang_index,"CDS")[:44])
                                                backbone_list[-1].append(backbone_list[-3][2])
                                                del backbone_list[-3]
                                            elif overhang in backbone_list[-3][0][bsa_index[x]+6:bsa_index[x]+37]:
                                                overhang_index = backbone_list[-3][0].find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                backbone_list[-2].append(silent_mut(backbone_list[-3][0][:overhang_index],bsa_index[x],"CDS")+overhang+"agagaccaa")
                                                backbone_list[-2].append(backbone_list[-3][1])
                                                backbone_list[-2].append(reverse_complement(silent_mut(backbone_list[-3][0][:overhang_index],bsa_index[x],"CDS")[-40:]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+backbone_list[-3][0][overhang_index:])
                                                backbone_list[-1].append("aaggtctca" + backbone_list[-3][0][overhang_index:overhang_index+44])
                                                backbone_list[-1].append(backbone_list[-3][2])
                                                del backbone_list[-3]
                                            else:
                                                pass
                            else:
                                backbone_list = [[new_backbone_sequence,new_backbone_sequence[:60],reverse_complement(new_backbone_sequence[-60:])]]
                            new_builds_list.append(new_unpacked_list)

                    #Perform scarless Golden Gate assembly computation
                    if golden_gate_method == "scarless_assembly":
                        new_backbone_sequence = app.registry.get(session_id)["new_backbone_sequence"]
                        new_builds_list = []
                        for unpacked_list in builds_list:
                            new_unpacked_list = []
                            contains_bsaI = False
                            for part in unpacked_list:
                                part.assembly_method = "Type_II_Restriction_Enzyme"
                                part.sequence = part.sequence.lower()
                                if "ggtctc" in part.sequence or "gagacc" in part.sequence or "ggtctc" in backbone_sequence or "gagacc" in backbone_sequence:
                                    golden_gate_error = "BsaI_in_seq" +"|"+ part.name
                                    contains_bsaI = True
                            #Given a list of linkers, try to find an efficient assignment of linkers for scarless GG assembly
                            def gg_scar_opt(unpacked_list,gg_overhangs):
                                new_backbone_sequence = app.registry.get(session_id)["new_backbone_sequence"]
                                golden_gate_overhangs = gg_overhangs
                                if remove_bsaI == "yes":
                                    golden_gate_error = ""
                                    bsaI_removal_combs = []
                                    bsaI_removal_key = []
                                    for part in unpacked_list:
                                        part.sequence = part.sequence.lower()
                                        bsa_index = []
                                        ii=0
                                        while len(bsa_index) != (part.sequence.count("ggtctc")+part.sequence.count("gagacc")) and ii<10*len(unpacked_list):
                                            if "ggtctc" in part.sequence and "gagacc" not in part.sequence:
                                                if bsa_index == []:
                                                    bsa_index.append(part.sequence.find("ggtctc"))
                                                else:
                                                    bsa_index.append(part.sequence.find("ggtctc",bsa_index[-1]+1))
                                            elif "gagacc" in part.sequence and "ggtctc" not in part.sequence:
                                                if bsa_index == []:
                                                    bsa_index.append(part.sequence.find("gagacc"))
                                                else:
                                                    bsa_index.append(part.sequence.find("gagacc",bsa_index[-1]+1))
                                            elif "ggtctc" in part.sequence and "gagacc" in part.sequence:
                                                if bsa_index == []:
                                                    bsa_index.append(min(part.sequence.find("ggtctc"),part.sequence.find("gagacc")))
                                                else:
                                                    if part.sequence.find("ggtctc",bsa_index[-1]+1) == -1 or part.sequence.find("gagacc",bsa_index[-1]+1) == -1:
                                                        bsa_index.append(max(part.sequence.find("ggtctc",bsa_index[-1]+1),part.sequence.find("gagacc",bsa_index[-1]+1)))
                                                    else:
                                                        bsa_index.append(min(part.sequence.find("ggtctc",bsa_index[-1]+1),part.sequence.find("gagacc",bsa_index[-1]+1)))
                                            ii+=1
                                        if "ggtctc" in part.sequence or "gagacc" in part.sequence:
                                            for xx in bsa_index:
                                                bsaI_removal_key.append(part)
                                            for i in range(len(bsa_index)):
                                                bsaI_removal_combs.append([])
                                                for overhang in golden_gate_overhangs:
                                                    if overhang in part.sequence[bsa_index[i]-30:bsa_index[i]+37] and ((part.sequence[bsa_index[i]-30:bsa_index[i]+37].find(overhang)+bsa_index[i]-30)<(bsa_index[i]-3) or (part.sequence[bsa_index[i]-30:bsa_index[i]+37].find(overhang)+bsa_index[i]-30)>=(bsa_index[i]+6)):
                                                        bsaI_removal_combs[-1].append(overhang)
                                    new_backbone_sequence = new_backbone_sequence.lower()
                                    bsa_index = []
                                    ii = 0
                                    while len(bsa_index) != (backbone_sequence.count("ggtctc")+backbone_sequence.count("gagacc")) and ii<10:
                                        if "ggtctc" in backbone_sequence and "gagacc" not in backbone_sequence:
                                            if bsa_index == []:
                                                bsa_index.append(backbone_sequence.find("ggtctc"))
                                            else:
                                                bsa_index.append(backbone_sequence.find("ggtctc",bsa_index[-1]+1))
                                        elif "gagacc" in backbone_sequence and "ggtctc" not in backbone_sequence:
                                            if bsa_index == []:
                                                bsa_index.append(backbone_sequence.find("gagacc"))
                                            else:
                                                bsa_index.append(backbone_sequence.find("gagacc",bsa_index[-1]+1))
                                        elif "ggtctc" in backbone_sequence and "gagacc" in backbone_sequence:
                                            if bsa_index == []:
                                                bsa_index.append(min(backbone_sequence.find("ggtctc"),backbone_sequence.find("gagacc")))
                                            else:
                                                if backbone_sequence.find("ggtctc",bsa_index[-1]+1) == -1 or backbone_sequence.find("gagacc",bsa_index[-1]+1) == -1:
                                                    bsa_index.append(max(backbone_sequence.find("ggtctc",bsa_index[-1]+1),backbone_sequence.find("gagacc",bsa_index[-1]+1)))
                                                else:
                                                    bsa_index.append(min(backbone_sequence.find("ggtctc",bsa_index[-1]+1),backbone_sequence.find("gagacc",bsa_index[-1]+1)))
                                        ii+=1
                                    if "ggtctc" in backbone_sequence or "gagacc" in backbone_sequence:
                                        for xx in bsa_index:
                                            bsaI_removal_key.append("backbone")
                                        for i in range(len(bsa_index)):
                                            bsaI_removal_combs.append([])
                                            for overhang in golden_gate_overhangs:
                                                if overhang in backbone_sequence[bsa_index[i]-30:bsa_index[i]+37] and ((backbone_sequence[bsa_index[i]-30:bsa_index[i]+37].find(overhang)+bsa_index[i]-30)<(bsa_index[i]-3) or (backbone_sequence[bsa_index[i]-30:bsa_index[i]+37].find(overhang)+bsa_index[i]-30)>=(bsa_index[i]+6)):
                                                    bsaI_removal_combs[-1].append(overhang)
                                else:
                                    bsaI_removal_combs = []
                                    bsaI_removal_key = []
                                return bsaI_removal_combs,bsaI_removal_key,golden_gate_error
                            #Find an efficient linker set for the scarless assembly, including removing internal BsaI sites
                            if remove_bsaI == "yes":
                                golden_gate_error = ""
                                breakit = False
                                for p in range(10,51):
                                    for q in range(5):
                                        with open('/var/www/ibiocad/iBioCAD/overhangsets/setsof%s.csv'%p,'r') as f:
                                            reader = csv.reader(f, delimiter=",")
                                            temp = list(reader)[1:]
                                        gg_overhangs = []
                                        for x in range(p):
                                            gg_overhangs.append(temp[q][x].lower())
                                        bsaI_removal_combs,bsaI_removal_key,golden_gate_error = gg_scar_opt(unpacked_list,gg_overhangs)
                                        gg_opt = golden_gate_optimization(unpacked_list,backbone_sequence,gg_overhangs,bsaI_combs=bsaI_removal_combs)
                                        if gg_opt is not None:
                                            breakit = True
                                        if breakit:
                                            break
                                    if breakit:
                                        break
                            #Find an efficient linker set for the scarless assembly if there are no internal BsaI sites
                            elif contains_bsaI==False:
                                golden_gate_error = ""
                                bsaI_removal_key = []
                                breakit = False
                                for p in range(10,51):
                                    for q in range(5):
                                        with open('/var/www/ibiocad/iBioCAD/overhangsets/setsof%s.csv'%p,'r') as f:
                                            reader = csv.reader(f, delimiter=",")
                                            temp = list(reader)[1:]
                                        gg_overhangs = []
                                        for x in range(p):
                                            gg_overhangs.append(temp[q][x].lower())
                                        gg_opt = golden_gate_optimization(unpacked_list,backbone_sequence,gg_overhangs)
                                        if gg_opt is not None:
                                            breakit = True
                                        if breakit:
                                            break
                                    if breakit:
                                        break
                            #Find an efficient linker set for scarless assembly, if internal BsaI sites should not be removed
                            else:
                                bsaI_removal_key = []
                            #Modify sequences and design primers
                            #First, modify sequences at the junctions between parts
                            for i in range(len(unpacked_list)):
                                if golden_gate_error != "" or remove_bsaI == "no":
                                    break
                                if bsaI_removal_key is not None:
                                    bsa_len = len(bsaI_removal_key)
                                else:
                                    bsa_len = 0
                                if i == 0:
                                    if gg_opt[0] in new_backbone_sequence[-35:]:
                                        new_backbone_sequence = new_backbone_sequence[:new_backbone_sequence.find(gg_opt[0],len(new_backbone_sequence)-35)] + gg_opt[0] + "agagaccaa"
                                        unpacked_list[0].primer_forward = "aaggtctca" + new_backbone_sequence[new_backbone_sequence.find(gg_opt[0],len(new_backbone_sequence)-35):-9] + unpacked_list[0].sequence[:40]
                                        unpacked_list[0].sequence = "aaggtctca" + new_backbone_sequence[new_backbone_sequence.find(gg_opt[0],len(new_backbone_sequence)-35):-9] + unpacked_list[0].sequence
                                    elif gg_opt[0] in unpacked_list[0].sequence[:35]:
                                        new_backbone_sequence = new_backbone_sequence + unpacked_list[0].sequence[:unpacked_list[0].sequence.find(gg_opt[0])] + gg_opt[0] + "agagaccaa"
                                        unpacked_list[0].primer_forward = "aaggtctca" + unpacked_list[0].sequence[unpacked_list[0].sequence.find(gg_opt[0]):unpacked_list[0].sequence.find(gg_opt[0])+40]
                                        unpacked_list[0].sequence = "aaggtctca" + unpacked_list[0].sequence[unpacked_list[0].sequence.find(gg_opt[0]):]
                                    if len(unpacked_list)==1:
                                        if gg_opt[-1-bsa_len] in new_backbone_sequence[:35]:
                                            unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[-1].sequence[-40:] + new_backbone_sequence[9:new_backbone_sequence.find(gg_opt[-1-bsa_len])] + gg_opt[-1-bsa_len] + "agagaccaa")
                                            unpacked_list[-1].sequence = unpacked_list[-1].sequence + new_backbone_sequence[:new_backbone_sequence.find(gg_opt[-1-bsa_len])] + gg_opt[-1-bsa_len] + "agagaccaa"
                                            new_backbone_sequence = "aaggtctca" + new_backbone_sequence[new_backbone_sequence.find(gg_opt[-1-bsa_len]):]
                                        elif gg_opt[-1-bsa_len] in unpacked_list[-1].sequence[-35:]:
                                            new_backbone_sequence = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-1-bsa_len],len(unpacked_list[-1].sequence)-35):] + new_backbone_sequence
                                            unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-1-bsa_len],len(unpacked_list[-1].sequence)-35)-40:unpacked_list[-1].sequence.find(gg_opt[-1-bsa_len],len(unpacked_list[-1].sequence)-35)] + gg_opt[-1-bsa_len] + "agagaccaa")
                                            unpacked_list[-1].sequence = unpacked_list[-1].sequence[:unpacked_list[-1].sequence.find(gg_opt[-1-bsa_len],len(unpacked_list[-1].sequence)-35)] + gg_opt[-1-bsa_len] + "agagaccaa"
                                        else:
                                            pass
                                        #if gg_opt[-2-bsa_len] in unpacked_list[-1].sequence[:35]:
                                        #    unpacked_list[-2].primer_reverse = reverse_complement(unpacked_list[-2].sequence[-40:] + unpacked_list[-1].sequence[:unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len])] + gg_opt[-2-bsa_len] + "agagaccaa")
                                        #    unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len]):unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len])+40]
                                        #    unpacked_list[-2].sequence = unpacked_list[-2].sequence + unpacked_list[-1].sequence[:unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len])] + gg_opt[-2-bsa_len] + "agagaccaa"
                                        #    unpacked_list[-1].sequence = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len]):]
                                        #elif gg_opt[-2-bsa_len] in unpacked_list[-2].sequence[-35:]:
                                        #    unpacked_list[-2].primer_reverse = reverse_complement(unpacked_list[-2].sequence[unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35)-40:unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35)] + gg_opt[-2-bsa_len] + "agagaccaa")
                                        #    unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[-2].sequence[unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35):] + unpacked_list[-1].sequence[:40]
                                        #    unpacked_list[-1].sequence = "aaggtctca" + unpacked_list[-2].sequence[unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35):] + unpacked_list[-1].sequence
                                        #    unpacked_list[-2].sequence = unpacked_list[-2].sequence[:unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35)] + gg_opt[-2-bsa_len] + "agagaccaa"
                                        #else:
                                        #    pass
                                elif i == len(unpacked_list)-1:
                                    if gg_opt[-1-bsa_len] in new_backbone_sequence[:35]:
                                        unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[-1].sequence[-40:] + new_backbone_sequence[9:new_backbone_sequence.find(gg_opt[-1-bsa_len])] + gg_opt[-1-bsa_len] + "agagaccaa")
                                        unpacked_list[-1].sequence = unpacked_list[-1].sequence + new_backbone_sequence[:new_backbone_sequence.find(gg_opt[-1-bsa_len])] + gg_opt[-1-bsa_len] + "agagaccaa"
                                        new_backbone_sequence = "aaggtctca" + new_backbone_sequence[new_backbone_sequence.find(gg_opt[-1-bsa_len]):]
                                    elif gg_opt[-1-bsa_len] in unpacked_list[-1].sequence[-35:]:
                                        new_backbone_sequence = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-1-bsa_len],len(unpacked_list[-1].sequence)-35):] + new_backbone_sequence
                                        unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-1-bsa_len],len(unpacked_list[-1].sequence)-35)-40:unpacked_list[-1].sequence.find(gg_opt[-1-bsa_len],len(unpacked_list[-1].sequence)-35)] + gg_opt[-1-bsa_len] + "agagaccaa")
                                        unpacked_list[-1].sequence = unpacked_list[-1].sequence[:unpacked_list[-1].sequence.find(gg_opt[-1-bsa_len],len(unpacked_list[-1].sequence)-35)] + gg_opt[-1-bsa_len] + "agagaccaa"
                                    else:
                                        pass
                                    if gg_opt[-2-bsa_len] in unpacked_list[-1].sequence[:35]:
                                        unpacked_list[-2].primer_reverse = reverse_complement(unpacked_list[-2].sequence[-40:] + unpacked_list[-1].sequence[:unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len])] + gg_opt[-2-bsa_len] + "agagaccaa")
                                        unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len]):unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len])+40]
                                        unpacked_list[-2].sequence = unpacked_list[-2].sequence + unpacked_list[-1].sequence[:unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len])] + gg_opt[-2-bsa_len] + "agagaccaa"
                                        unpacked_list[-1].sequence = "aaggtctca" + unpacked_list[-1].sequence[unpacked_list[-1].sequence.find(gg_opt[-2-bsa_len]):]
                                    elif gg_opt[-2-bsa_len] in unpacked_list[-2].sequence[-35:]:
                                        unpacked_list[-2].primer_reverse = reverse_complement(unpacked_list[-2].sequence[unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35)-40:unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35)] + gg_opt[-2-bsa_len] + "agagaccaa")
                                        unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[-2].sequence[unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35):] + unpacked_list[-1].sequence[:40]
                                        unpacked_list[-1].sequence = "aaggtctca" + unpacked_list[-2].sequence[unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35):] + unpacked_list[-1].sequence
                                        unpacked_list[-2].sequence = unpacked_list[-2].sequence[:unpacked_list[-2].sequence.find(gg_opt[-2-bsa_len],len(unpacked_list[-2].sequence)-35)] + gg_opt[-2-bsa_len] + "agagaccaa"
                                    else:
                                        pass
                                else:
                                    if gg_opt[i] in unpacked_list[i].sequence[:35]:
                                        unpacked_list[i-1].primer_reverse = reverse_complement(unpacked_list[i-1].sequence[-40:] + unpacked_list[i].sequence[:unpacked_list[i].sequence.find(gg_opt[i])] + gg_opt[i] + "agagaccaa")
                                        unpacked_list[i].primer_forward = "aaggtctca" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(gg_opt[i]):unpacked_list[i].sequence.find(gg_opt[i])+40]
                                        unpacked_list[i-1].sequence = unpacked_list[i-1].sequence + unpacked_list[i].sequence[:unpacked_list[i].sequence.find(gg_opt[i])] + gg_opt[i] + "agagaccaa"
                                        unpacked_list[i].sequence = "aaggtctca" + unpacked_list[i].sequence[unpacked_list[i].sequence.find(gg_opt[i]):]
                                    elif gg_opt[i] in unpacked_list[i-1].sequence[-35:]:
                                        unpacked_list[i-1].primer_reverse = reverse_complement(unpacked_list[i-1].sequence[unpacked_list[i-1].sequence.find(gg_opt[i],len(unpacked_list[i-1].sequence)-35)-40:unpacked_list[i-1].sequence.find(gg_opt[i],len(unpacked_list[i-1].sequence)-35)] + gg_opt[i] + "agagaccaa")
                                        unpacked_list[i].primer_forward = "aaggtctca" + unpacked_list[i-1].sequence[unpacked_list[i-1].sequence.find(gg_opt[i],len(unpacked_list[i-1].sequence)-35):] + unpacked_list[i].sequence[:40]
                                        unpacked_list[i].sequence = "aaggtctca" + unpacked_list[i-1].sequence[unpacked_list[i-1].sequence.find(gg_opt[i],len(unpacked_list[i-1].sequence)-35):] + unpacked_list[i].sequence
                                        unpacked_list[i-1].sequence = unpacked_list[i-1].sequence[:unpacked_list[i-1].sequence.find(gg_opt[i],len(unpacked_list[i-1].sequence)-35)] + gg_opt[i] + "agagaccaa"
                                    else:
                                        pass
                            #Then, remove internal BsaI sites
                            for i in range(len(unpacked_list)):
                                if golden_gate_error != "" or remove_bsaI == "no":
                                    break
                                if unpacked_list[i] in bsaI_removal_key:
                                    bsa_index = []
                                    ii=0
                                    while len(bsa_index) != (unpacked_list[i].sequence.count("ggtctc")+unpacked_list[i].sequence.count("gagacc"))-2 and ii<10*len(unpacked_list):
                                        if "ggtctc" in unpacked_list[i].sequence[13:] and "gagacc" not in unpacked_list[i].sequence[:len(unpacked_list[i].sequence)-13]:
                                            if bsa_index == []:
                                                bsa_index.append(unpacked_list[i].sequence.find("ggtctc",13))
                                            else:
                                                bsa_index.append(unpacked_list[i].sequence.find("ggtctc",bsa_index[-1]+1))
                                        elif "gagacc" in unpacked_list[i].sequence[:len(unpacked_list[i].sequence)-13] and "ggtctc" not in unpacked_list[i].sequence[13:]:
                                            if bsa_index == []:
                                                bsa_index.append(unpacked_list[i].sequence.find("gagacc",0,len(unpacked_list[i].sequence)-13))
                                            else:
                                                bsa_index.append(unpacked_list[i].sequence.find("gagacc",bsa_index[-1]+1))
                                        elif "ggtctc" in unpacked_list[i].sequence[13:] and "gagacc" in unpacked_list[i].sequence[:len(unpacked_list[i].sequence)-13]:
                                            if bsa_index == []:
                                                bsa_index.append(min(unpacked_list[i].sequence.find("ggtctc",13),unpacked_list[i].sequence.find("gagacc",0,len(unpacked_list[i].sequence)-13)))
                                            else:
                                                bsa_index.append(min(unpacked_list[i].sequence.find("ggtctc",bsa_index[-1]+1),unpacked_list[i].sequence.find("gagacc",bsa_index[-1]+1)))
                                        ii+=1
                                    for x in range(bsaI_removal_key.count(unpacked_list[i])):
                                        part_index = bsaI_removal_key.index(unpacked_list[i])+x
                                        overhang = gg_opt[len(unpacked_list)+1+part_index]
                                        if x == 0:
                                            if bsa_index[x]-30 < 0:
                                                if overhang in unpacked_list[i].sequence[0:bsa_index[x]]:
                                                    overhang_index = unpacked_list[i].sequence.find(overhang,0,bsa_index[x])
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"a",unpacked_list[i].type,unpacked_list[i].sequence[:overhang_index]+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = unpacked_list[i].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[i].sequence[overhang_index-40:overhang_index]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"b",unpacked_list[i].type,"aaggtctca"+silent_mut(unpacked_list[i].sequence[overhang_index:],bsa_index[x]-overhang_index,unpacked_list[i].type)))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + silent_mut(unpacked_list[i].sequence[overhang_index:],bsa_index[x]-overhang_index,unpacked_list[i].type)[:44]
                                                    new_unpacked_list[-1].primer_reverse = unpacked_list[i].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                elif overhang in unpacked_list[i].sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                    overhang_index = unpacked_list[i].sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"a",unpacked_list[i].type,silent_mut(unpacked_list[i].sequence[:overhang_index],bsa_index[x],unpacked_list[i].type)+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = unpacked_list[i].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(silent_mut(unpacked_list[i].sequence[:overhang_index],bsa_index[x],unpacked_list[i].type)[-40:]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"b",unpacked_list[i].type,"aaggtctca"+unpacked_list[i].sequence[overhang_index:]))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[i].sequence[overhang_index:overhang_index+44]
                                                    new_unpacked_list[-1].primer_reverse = unpacked_list[i].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                else:
                                                    pass
                                            else:
                                                if overhang in unpacked_list[i].sequence[bsa_index[x]-30:bsa_index[x]]:
                                                    overhang_index = unpacked_list[i].sequence.find(overhang,bsa_index[x]-30,bsa_index[x])
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"a",unpacked_list[i].type,unpacked_list[i].sequence[:overhang_index]+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = unpacked_list[i].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(unpacked_list[i].sequence[overhang_index-40:overhang_index]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"b",unpacked_list[i].type,"aaggtctca"+silent_mut(unpacked_list[i].sequence[overhang_index:],bsa_index[x]-overhang_index,unpacked_list[i].type)))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + silent_mut(unpacked_list[i].sequence[overhang_index:],bsa_index[x]-overhang_index,unpacked_list[i].type)[:44]
                                                    new_unpacked_list[-1].primer_reverse = unpacked_list[i].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                elif overhang in unpacked_list[i].sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                    overhang_index = unpacked_list[i].sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"a",unpacked_list[i].type,silent_mut(unpacked_list[i].sequence[:overhang_index],bsa_index[x],unpacked_list[i].type)+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = unpacked_list[i].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(silent_mut(unpacked_list[i].sequence[:overhang_index],bsa_index[x],unpacked_list[i].type)[-40:]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(unpacked_list[i].name+"-"+str(x+1)+"b",unpacked_list[i].type,"aaggtctca"+unpacked_list[i].sequence[overhang_index:]))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + unpacked_list[i].sequence[overhang_index:overhang_index+44]
                                                    new_unpacked_list[-1].primer_reverse = unpacked_list[i].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                else:
                                                    pass
                                        else:
                                            if bsa_index[x]-30 < 0:
                                                if overhang in new_unpacked_list[-1].sequence[0:bsa_index[x]]:
                                                    overhang_index = new_unpacked_list[-1].sequence.find(overhang,0,bsa_index[x])
                                                    new_unpacked_list.append(Part(new_unpacked_list[-1].name+"-"+str(x+1)+"a",new_unpacked_list[-1].type,new_unpacked_list[-1].sequence[:overhang_index]+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = new_unpacked_list[-2].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(new_unpacked_list[-2].sequence[overhang_index-40:overhang_index]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(new_unpacked_list[-2].name+"-"+str(x+1)+"b",new_unpacked_list[-2].type,"aaggtctca"+silent_mut(new_unpacked_list[-2].sequence[overhang_index:],bsa_index[x]-overhang_index,new_unpacked_list[-2].type)))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + silent_mut(new_unpacked_list[-3].sequence[overhang_index:],bsa_index[x]-overhang_index,new_unpacked_list[-3].type)[:44]
                                                    new_unpacked_list[-1].primer_reverse = new_unpacked_list[-3].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    del new_unpacked_list[-3]
                                                elif overhang in new_unpacked_list[-1].sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                    overhang_index = new_unpacked_list[-1].sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                    new_unpacked_list.append(Part(new_unpacked_list[-1].name+"-"+str(x+1)+"a",new_unpacked_list[-1].type,silent_mut(new_unpacked_list[-1].sequence[:overhang_index],bsa_index[x],new_unpacked_list[-1].type)+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = new_unpacked_list[-2].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(silent_mut(new_unpacked_list[-2].sequence[:overhang_index],bsa_index[x],new_unpacked_list[-2].type)[-40:]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(new_unpacked_list[-2].name+"-"+str(x+1)+"b",new_unpacked_list[-2].type,"aaggtctca"+new_unpacked_list[-2].sequence[overhang_index:]))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + new_unpacked_list[-3].sequence[overhang_index:overhang_index+44]
                                                    new_unpacked_list[-1].primer_reverse = new_unpacked_list[-3].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    del new_unpacked_list[-3]
                                                else:
                                                    pass
                                            else:
                                                if overhang in new_unpacked_list[-1].sequence[bsa_index[x]-30:bsa_index[x]]:
                                                    overhang_index = new_unpacked_list[-1].sequence.find(overhang,bsa_index[x]-30,bsa_index[x])
                                                    new_unpacked_list.append(Part(new_unpacked_list[-1].name+"-"+str(x+1)+"a",new_unpacked_list[-1].type,new_unpacked_list[-1].sequence[:overhang_index]+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = new_unpacked_list[-2].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(new_unpacked_list[-2].sequence[overhang_index-40:overhang_index]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(new_unpacked_list[-2].name+"-"+str(x+1)+"b",new_unpacked_list[-2].type,"aaggtctca"+silent_mut(new_unpacked_list[-2].sequence[overhang_index:],bsa_index[x]-overhang_index,new_unpacked_list[-2].type)))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + silent_mut(new_unpacked_list[-3].sequence[overhang_index:],bsa_index[x]-overhang_index,new_unpacked_list[-3].type)[:44]
                                                    new_unpacked_list[-1].primer_reverse = new_unpacked_list[-3].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    del new_unpacked_list[-3]
                                                elif overhang in new_unpacked_list[-1].sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                    overhang_index = new_unpacked_list[-1].sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                    new_unpacked_list.append(Part(new_unpacked_list[-1].name+"-"+str(x+1)+"a",new_unpacked_list[-1].type,silent_mut(new_unpacked_list[-1].sequence[:overhang_index],bsa_index[x],new_unpacked_list[-1].type)+overhang+"agagaccaa"))
                                                    new_unpacked_list[-1].primer_forward = new_unpacked_list[-2].primer_forward
                                                    new_unpacked_list[-1].primer_reverse = reverse_complement(silent_mut(new_unpacked_list[-2].sequence[:overhang_index],bsa_index[x],new_unpacked_list[-2].type)[-40:]+overhang+"agagaccaa")
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    if len(bsa_index)>x+1:
                                                        for q in range(x+1,len(bsa_index)):
                                                            bsa_index[q] = bsa_index[q]-len(new_unpacked_list[-1].sequence[:-9])+13
                                                    new_unpacked_list.append(Part(new_unpacked_list[-2].name+"-"+str(x+1)+"b",new_unpacked_list[-2].type,"aaggtctca"+new_unpacked_list[-2].sequence[overhang_index:]))
                                                    new_unpacked_list[-1].primer_forward = "aaggtctca" + new_unpacked_list[-3].sequence[overhang_index:overhang_index+44]
                                                    new_unpacked_list[-1].primer_reverse = new_unpacked_list[-3].primer_reverse
                                                    new_unpacked_list[-1].assembly_method = "Type_II_Restriction_Enzyme"
                                                    del new_unpacked_list[-3]
                                                else:
                                                    pass
                                else:
                                    new_unpacked_list.append(unpacked_list[i])
                            #Remove BsaI sites internal in backbone sequence
                            if "backbone" in bsaI_removal_key:
                                bsa_index = []
                                ii=0
                                while len(bsa_index) != (new_backbone_sequence.count("ggtctc")+new_backbone_sequence.count("gagacc"))-2 and ii<10:
                                    if "ggtctc" in new_backbone_sequence[13:] and "gagacc" not in new_backbone_sequence[:len(new_backbone_sequence)-13]:
                                        if bsa_index == []:
                                            bsa_index.append(new_backbone_sequence.find("ggtctc",13))
                                        else:
                                            bsa_index.append(new_backbone_sequence.find("ggtctc",bsa_index[-1]+1))
                                    elif "gagacc" in new_backbone_sequence[:len(new_backbone_sequence)-13] and "ggtctc" not in new_backbone_sequence[13:]:
                                        if bsa_index == []:
                                            bsa_index.append(new_backbone_sequence.find("gagacc",0,len(new_backbone_sequence)-13))
                                        else:
                                            bsa_index.append(new_backbone_sequence.find("gagacc",bsa_index[-1]+1))
                                    elif "ggtctc" in new_backbone_sequence[13:] and "gagacc" in new_backbone_sequence[:len(new_backbone_sequence)-13]:
                                        if bsa_index == []:
                                            bsa_index.append(min(new_backbone_sequence.find("ggtctc",13),new_backbone_sequence.find("gagacc",0,len(new_backbone_sequence)-13)))
                                        else:
                                            bsa_index.append(min(new_backbone_sequence.find("ggtctc",bsa_index[-1]+1),new_backbone_sequence.find("gagacc",bsa_index[-1]+1)))
                                    ii+=1
                                backbone_list = []
                                new_backbone_primers = [new_backbone_sequence[:40],reverse_complement(new_backbone_sequence[-40:])]
                                for x in range(bsaI_removal_key.count("backbone")):
                                    backbone_list.append([])
                                    backbone_list.append([])
                                    part_index = bsaI_removal_key.index("backbone")+x
                                    overhang = gg_opt[len(unpacked_list)+1+part_index]
                                    if x == 0:
                                        if bsa_index[x]-30 < 0:
                                            if overhang in new_backbone_sequence[0:bsa_index[x]]:
                                                overhang_index = new_backbone_sequence.find(overhang,0,bsa_index[x])
                                                backbone_list[-2].append(new_backbone_sequence[:overhang_index]+overhang+"agagaccaa")
                                                backbone_list[-2].append(new_backbone_primers[0])
                                                backbone_list[-2].append(reverse_complement(new_backbone_sequence[overhang_index-40:overhang_index]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+silent_mut(new_backbone_sequence[overhang_index:],bsa_index[x]-overhang_index,"CDS"))
                                                backbone_list[-1].append("aaggtctca" + silent_mut(new_backbone_sequence[overhang_index:],bsa_index[x]-overhang_index,"CDS")[:44])
                                                backbone_list[-1].append(new_backbone_primers[1])
                                            elif overhang in new_backbone_sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                overhang_index = new_backbone_sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                backbone_list[-2].append(silent_mut(new_backbone_sequence[:overhang_index],bsa_index[x],"CDS")+overhang+"agagaccaa")
                                                backbone_list[-2].append(new_backbone_primers[0])
                                                backbone_list[-2].append(reverse_complement(silent_mut(new_backbone_sequence[:overhang_index],bsa_index[x],"CDS")[-40:]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+new_backbone_sequence[overhang_index:])
                                                backbone_list[-1].append("aaggtctca" + new_backbone_sequence[overhang_index:overhang_index+44])
                                                backbone_list[-1].append(new_backbone_primers[1])
                                            else:
                                                pass
                                        else:
                                            if overhang in new_backbone_sequence[bsa_index[x]-30:bsa_index[x]]:
                                                overhang_index = new_backbone_sequence.find(overhang,bsa_index[x]-30,bsa_index[x])
                                                backbone_list[-2].append(new_backbone_sequence[:overhang_index]+overhang+"agagaccaa")
                                                backbone_list[-2].append(new_backbone_primers[0])
                                                backbone_list[-2].append(reverse_complement(new_backbone_sequence[overhang_index-40:overhang_index]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+silent_mut(new_backbone_sequence[overhang_index:],bsa_index[x]-overhang_index,"CDS"))
                                                backbone_list[-1].append("aaggtctca" + silent_mut(new_backbone_sequence[overhang_index:],bsa_index[x]-overhang_index,"CDS")[:44])
                                                backbone_list[-1].append(new_backbone_primers[1])
                                            elif overhang in new_backbone_sequence[bsa_index[x]+6:bsa_index[x]+37]:
                                                overhang_index = new_backbone_sequence.find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                backbone_list[-2].append(silent_mut(new_backbone_sequence[:overhang_index],bsa_index[x],"CDS")+overhang+"agagaccaa")
                                                backbone_list[-2].append(new_backbone_primers[0])
                                                backbone_list[-2].append(reverse_complement(silent_mut(new_backbone_sequence[:overhang_index],bsa_index[x],"CDS")[-40:]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+new_backbone_sequence[overhang_index:])
                                                backbone_list[-1].append("aaggtctca" + new_backbone_sequence[overhang_index:overhang_index+44])
                                                backbone_list[-1].append(new_backbone_primers[1])
                                            else:
                                                pass
                                    else:
                                        if bsa_index[x]-30 < 0:
                                            if overhang in backbone_list[-3][0][0:bsa_index[x]]:
                                                overhang_index = backbone_list[-3][0].find(overhang,0,bsa_index[x])
                                                backbone_list[-2].append(backbone_list[-3][0][:overhang_index]+overhang+"agagaccaa")
                                                backbone_list[-2].append(backbone_list[-3][1])
                                                backbone_list[-2].append(reverse_complement(backbone_list[-3][0][overhang_index-40:overhang_index]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+silent_mut(backbone_list[-3][0][overhang_index:],bsa_index[x]-overhang_index,"CDS"))
                                                backbone_list[-1].append("aaggtctca" + silent_mut(backbone_list[-3][0][overhang_index:],bsa_index[x]-overhang_index,"CDS")[:44])
                                                backbone_list[-1].append(backbone_list[-3][2])
                                                del backbone_list[-3]
                                            elif overhang in backbone_list[-3][0][bsa_index[x]+6:bsa_index[x]+37]:
                                                overhang_index = backbone_list[-3][0].find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                backbone_list[-2].append(silent_mut(backbone_list[-3][0][:overhang_index],bsa_index[x],"CDS")+overhang+"agagaccaa")
                                                backbone_list[-2].append(backbone_list[-3][1])
                                                backbone_list[-2].append(reverse_complement(silent_mut(backbone_list[-3][0][:overhang_index],bsa_index[x],"CDS")[-40:]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+backbone_list[-3][0][overhang_index:])
                                                backbone_list[-1].append("aaggtctca" + backbone_list[-3][0][overhang_index:overhang_index+44])
                                                backbone_list[-1].append(backbone_list[-3][2])
                                                del backbone_list[-3]
                                            else:
                                                pass
                                        else:
                                            if overhang in backbone_list[-3][0][bsa_index[x]-30:bsa_index[x]]:
                                                overhang_index = backbone_list[-3][0].find(overhang,bsa_index[x]-30,bsa_index[x])
                                                backbone_list[-2].append(backbone_list[-3][0][:overhang_index]+overhang+"agagaccaa")
                                                backbone_list[-2].append(backbone_list[-3][1])
                                                backbone_list[-2].append(reverse_complement(backbone_list[-3][0][overhang_index-40:overhang_index]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+silent_mut(backbone_list[-3][0][overhang_index:],bsa_index[x]-overhang_index,"CDS"))
                                                backbone_list[-1].append("aaggtctca" + silent_mut(backbone_list[-3][0][overhang_index:],bsa_index[x]-overhang_index,"CDS")[:44])
                                                backbone_list[-1].append(backbone_list[-3][2])
                                                del backbone_list[-3]
                                            elif overhang in backbone_list[-3][0][bsa_index[x]+6:bsa_index[x]+37]:
                                                overhang_index = backbone_list[-3][0].find(overhang,bsa_index[x]+6,bsa_index[x]+37)
                                                backbone_list[-2].append(silent_mut(backbone_list[-3][0][:overhang_index],bsa_index[x],"CDS")+overhang+"agagaccaa")
                                                backbone_list[-2].append(backbone_list[-3][1])
                                                backbone_list[-2].append(reverse_complement(silent_mut(backbone_list[-3][0][:overhang_index],bsa_index[x],"CDS")[-40:]+overhang+"agagaccaa"))
                                                if len(bsa_index)>x+1:
                                                    for q in range(x+1,len(bsa_index)):
                                                        bsa_index[q] = bsa_index[q]-len(backbone_list[-2][0][:-9])+13
                                                backbone_list[-1].append("aaggtctca"+backbone_list[-3][0][overhang_index:])
                                                backbone_list[-1].append("aaggtctca" + backbone_list[-3][0][overhang_index:overhang_index+44])
                                                backbone_list[-1].append(backbone_list[-3][2])
                                                del backbone_list[-3]
                                            else:
                                                pass
                            else:
                                backbone_list = [[new_backbone_sequence,new_backbone_sequence[:60],reverse_complement(new_backbone_sequence[-60:])]]
                            new_builds_list.append(new_unpacked_list)

                    #Optimize primer Tm based off config settings
                    if primer_optimization == "range":
                        for unpacked_list in new_builds_list:
                            for part in unpacked_list:
                                if mt.Tm_NN(part.primer_forward) < primer_tm_range[0]:
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                else:
                                    for i in range(len(part.primer_forward),32,-1):
                                        if mt.Tm_NN(part.primer_forward[:i]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i]) <= primer_tm_range[1]:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                            part.primer_forward = part.primer_forward[:i]
                                            break
                                        else:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                if mt.Tm_NN(part.primer_reverse) < primer_tm_range[0]:
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                else:
                                    for i in range(len(part.primer_reverse),32,-1):
                                        if mt.Tm_NN(part.primer_reverse[:i]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:i]) <= primer_tm_range[1]:
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:i])
                                            part.primer_reverse = part.primer_reverse[:i]
                                            break
                                        else:
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        breakit=False
                        for x in range(60,32,-1):
                            for z in range(60,32,-1):
                                if mt.Tm_NN(new_backbone_sequence[:x])>=primer_tm_range[0] and mt.Tm_NN(new_backbone_sequence[:x])<=primer_tm_range[1] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:])[:z]) >= primer_tm_range[0] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:])[:z]) <= primer_tm_range[1]:
                                    backbone_primers = [new_backbone_sequence[:x],reverse_complement(new_backbone_sequence[-46:])[:z]]
                                    backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:x]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:])[:z])]
                                    breakit=True
                                if breakit:
                                    break
                            if breakit:
                                break
                        if not breakit:
                            backbone_primers = [new_backbone_sequence[:60],reverse_complement(new_backbone_sequence[-60:])]
                            backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:60]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:]))]
                    elif primer_optimization == "near":
                        for unpacked_list in new_builds_list:
                            for part in unpacked_list:
                                breakit=False
                                for i in range(len(part.primer_forward),32,-1):
                                    for j in range(len(part.primer_reverse),32,-1):
                                        if abs(mt.Tm_NN(part.primer_forward[:i])-mt.Tm_NN(part.primer_reverse[:j])) <=5:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                            part.primer_forward = part.primer_forward[:i]
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j])
                                            part.primer_reverse = part.primer_reverse[:j]
                                            breakit=True
                                        if breakit:
                                            break
                                    if breakit:
                                        break
                                if not breakit:
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        breakit=False
                        for x in range(60,32,-1):
                            for z in range(60,32,-1):
                                if abs(mt.Tm_NN(new_backbone_sequence[:x])-mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:])[:z]))<=5:
                                    backbone_primers = [new_backbone_sequence[:x],reverse_complement(new_backbone_sequence[-60:])[:z]]
                                    backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:x]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:])[:z])]
                                    breakit=True
                                if breakit:
                                    break
                            if breakit:
                                break
                        if not breakit:
                            backbone_primers = [new_backbone_sequence[:60],reverse_complement(new_backbone_sequence[-60:])]
                            backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:60]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:]))]
                    elif primer_optimization == "both":
                        for unpacked_list in new_builds_list:
                            for part in unpacked_list:
                                if mt.Tm_NN(part.primer_forward) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse) < primer_tm_range[0]:
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                elif mt.Tm_NN(part.primer_forward) < primer_tm_range[0] and mt.Tm_NN(part.primer_reverse) >= primer_tm_range[0]:
                                    breakit=False
                                    part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                    for j in range(len(part.primer_reverse),32,-1):
                                        if abs(part.primer_forward_tm-mt.Tm_NN(part.primer_reverse[:j])) <=5:
                                            part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j])
                                            part.primer_reverse = part.primer_reverse[:j]
                                            breakit=True
                                    if not breakit:
                                        part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                elif mt.Tm_NN(part.primer_reverse) < primer_tm_range[0] and mt.Tm_NN(part.primer_forward) >= primer_tm_range[0]:
                                    breakit=False
                                    part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                                    for i in range(len(part.primer_forward),32,-1):
                                        if abs(part.primer_reverse_tm-mt.Tm_NN(part.primer_forward[:i])) <=5:
                                            part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                            part.primer_forward = part.primer_forward[:i]
                                            breakit=True
                                    if not breakit:
                                        part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                else:
                                    breakit=False
                                    for i in range(len(part.primer_forward),32,-1):
                                        for j in range(len(part.primer_reverse),32,-1):
                                            if abs(mt.Tm_NN(part.primer_forward[:i])-mt.Tm_NN(part.primer_reverse[:j]))<=5 and mt.Tm_NN(part.primer_forward[:i]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_forward[:i]) <= primer_tm_range[1] and mt.Tm_NN(part.primer_reverse[:j]) >= primer_tm_range[0] and mt.Tm_NN(part.primer_reverse[:j]) <= primer_tm_range[1]:
                                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward[:i])
                                                part.primer_forward = part.primer_forward[:i]
                                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse[:j])
                                                part.primer_reverse = part.primer_reverse[:j]
                                                breakit=True
                                            if breakit:
                                                break
                                        if breakit:
                                            break
                                    if not breakit:
                                        part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                        part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        breakit=False
                        for x in range(60,32,-1):
                            for z in range(60,32,-1):
                                if abs(mt.Tm_NN(new_backbone_sequence[:x])-mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:])[:z]))<=5 and mt.Tm_NN(new_backbone_sequence[:x])>=primer_tm_range[0] and mt.Tm_NN(new_backbone_sequence[:x])<=primer_tm_range[1] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:])[:z]) >= primer_tm_range[0] and mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:])[:z]) <= primer_tm_range[1]:
                                    backbone_primers = [new_backbone_sequence[:x],reverse_complement(new_backbone_sequence[-60:])[:z]]
                                    backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:x]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:])[:z])]
                                    breakit=True
                                if breakit:
                                    break
                            if breakit:
                                break
                        if not breakit:
                            backbone_primers = [new_backbone_sequence[:60],reverse_complement(new_backbone_sequence[-60:])]
                            backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:60]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:]))]
                    else:
                        for unpacked_list in new_builds_list:
                            for part in unpacked_list:
                                part.primer_forward_tm = mt.Tm_NN(part.primer_forward)
                                part.primer_reverse_tm = mt.Tm_NN(part.primer_reverse)
                        backbone_primers = [new_backbone_sequence[:60],reverse_complement(new_backbone_sequence[-60:])]
                        backbone_primers_tm = [mt.Tm_NN(new_backbone_sequence[:60]),mt.Tm_NN(reverse_complement(new_backbone_sequence[-60:]))]

                    #Register new parts and direct to output page
                    app.registry[session_id]["backbone_primers"] = backbone_primers
                    app.registry[session_id]["backbone_primers_tm"] = backbone_primers_tm
                    app.registry[session_id]["backbone_list"] = backbone_list
                    parts_list,session_id = self.update_part_list(updated_parts_list=parts_list)
                    app.registry[session_id]['builds_list'] = new_builds_list
                    app.registry[session_id]["new_backbone_sequence"] = new_backbone_sequence
                    if golden_gate_error == "":
                        self.redirect("/assembly")
                self.render("main_page.html",golden_gate_error=golden_gate_error,css=css,js=js,parts_list=parts_list)
            except Exception as e:
                logging.error(e)
                self.redirect("/error")


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
            try:
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
            except Exception as e:
                logging.error(e)
                self.redirect("/error")
        def post(self):
            try:
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
                                nameerror=""
                                for part in parts_list:
                                    if part.name == record.name:
                                        nameerror="part_already_input"
                                if nameerror=="":
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
            except Exception as e:
                logging.error(e)
                self.redirect("/error")

    #Write the results to CSV and to the assembly output page
    class AssemblyHandler(Handler):
        def get(self):
            try:
	            parts_list,session_id = self.get_parts_list()
	            application = webapp2.get_app()
	            builds_list = application.registry.get(session_id)['builds_list']
	            new_backbone_sequence = application.registry.get(session_id)["new_backbone_sequence"]
	            backbone_primers = application.registry.get(session_id)["backbone_primers"]
	            backbone_primers_tm = application.registry.get(session_id)["backbone_primers_tm"]
	            backbone_list = application.registry.get(session_id)["backbone_list"]
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
	                        csvdictwriter2 = csv.DictWriter(csvfile,fieldnames=['Name','','Sequence','Tm'])
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
	                    for num,backbone_part in enumerate(backbone_list,1):
	                        csvwriter.writerow(['Backbone'+str(num),'',backbone_part[0],''])
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
	                            csvwriter.writerow(['Name','','Sequence','Tm'])
	                        elif unpacked_list[0].assembly_method == "LCR":
	                            csvwriter.writerow(['Bridge','','Sequence'])
	                    if unpacked_list[0].assembly_method == "Yeast_Assembly" or unpacked_list[0].assembly_method == "Gibson_Assembly" or unpacked_list[0].assembly_method == "Type_II_Restriction_Enzyme":
	                        for num,backbone_part in enumerate(backbone_list,1):
	                            csvdictwriter2.writerow({'Name':("Backbone,forward"+str(num)),'':'','Sequence':backbone_part[1],'Tm':mt.Tm_NN(backbone_part[1])})
	                            csvdictwriter2.writerow({'Name':("Backbone,reverse"+str(num)),'':'','Sequence':backbone_part[2],'Tm':mt.Tm_NN(backbone_part[2])})
	                        for part in unpacked_list:
	                            csvdictwriter2.writerow({'Name':(part.name+",forward"),'':'','Sequence':part.primer_forward,'Tm':part.primer_forward_tm})
	                            csvdictwriter2.writerow({'Name':(part.name+",reverse"),'':'','Sequence':part.primer_reverse,'Tm':part.primer_reverse_tm})
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
	            self.render("assembly_page.html",builds_list=builds_list,part_dir=part_dir,new_backbone_sequence=new_backbone_sequence,backbone_primers=backbone_primers,backbone_primers_tm=backbone_primers_tm,backbone_list=backbone_list,data_uri=data_uri,filename=filename,css=css)
            except Exception as e:
               logging.error(e)
               self.redirect("/error")

    #Download XML file
    class ConstructDownloadHandler(Handler):
        def get(self):
            try:
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
            except Exception as e:
                logging.error(e)
                self.redirect("/error")

    #Read in config settings from config page
    class ConfigHandler(Handler):
        def get(self):
            try:
                if is_local:
                    #path on local
                    for record in SeqIO.parse("templates/pET-26b.fa","fasta"):
                        default_backbone = record
                else:
                    #path on server
                    for record in SeqIO.parse("/var/www/ibiocad/iBioCAD/templates/pET-26b.fa","fasta"):
                        default_backbone = record
                default_config = {"backbone":default_backbone,"Golden_gate_method":"scarless_assembly","backbone_primers_tm":[mt.Tm_NN(default_backbone.seq[:40]),mt.Tm_NN(default_backbone.seq[-40:])],"backbone_primers":[default_backbone.seq[:40],default_backbone.seq[-40:]]
                                  ,"Primer_optimization":"both","primer_tm":[52,60]}
                parts_list,session_id = self.get_parts_list()
                application = webapp2.get_app()
                if "assembly_config" not in application.registry.get(session_id).keys() or application.registry[session_id]["assembly_config"] is None:
                    application.registry[session_id]["assembly_config"] = default_config
                assembly_config = application.registry.get(session_id)["assembly_config"]
                self.render("config.html",css=css,assembly_config=assembly_config)
            except Exception as e:
                logging.error(e)
                self.redirect("/error")
        def post(self):
            try:
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
                        elif key=="primer_tm":
                            application.registry[session_id]["assembly_config"][key] = [int(self.request.POST.get("primer_tm0")),int(self.request.POST.get("primer_tm1"))]
                        else:
                            application.registry[session_id]["assembly_config"][key] = self.request.POST.get(key)
                    self.redirect("/")
            except Exception as e:
                logging.error(e)
                self.redirect("/error")

    #Render error page
    class ErrorHandler(Handler):
        def get(self):
            self.render("error_page.html",css=css)
        def post(self):
            self.redirect("/")

    #Initialize the web framework. URL paths are attached to their respective handlers
    #WSGI app must have the variable name "application" I think
    application = webapp2.WSGIApplication([
        ('/', MainHandler),('/break',breakHandler),('/inputpart',InputPartHandler),('/assembly',AssemblyHandler),('/construct_download',ConstructDownloadHandler),('/config',ConfigHandler),('/error',ErrorHandler)#,('/ExamplePath',exampleHandler)
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
except Exception as e:
    logging.error(e)
    pass