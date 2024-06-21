## Cargamos los paquetes
#import openpyxl
import requests
import cobra
import os
from cobra.io import read_sbml_model
from cobra.medium import minimal_medium
import warnings
warnings.filterwarnings("ignore")
#import json
#For FLYCOP
import smac
import logging
import numpy as np
from ConfigSpace.hyperparameters import UniformFloatHyperparameter
from ConfigSpace.hyperparameters import CategoricalHyperparameter #para la lista de nombres de modelos
from ConfigSpace.forbidden import ForbiddenEqualsClause, ForbiddenAndConjunction #para que se escojan 3 organismos diferentes en cada iteración de FLYCOP
from ConfigSpace.hyperparameters import Constant
# Import ConfigSpace and different types of parameters
from ConfigSpace import ConfigurationSpace, Categorical, Float, Configuration
from smac import HyperparameterOptimizationFacade
# Import SMAC-utilities
from smac import Scenario
import pandas as pd
import time
import cometspy as c
import os
import shutil
import cobra
from cobra import Reaction,Metabolite
import requests
from os.path import exists
import cobra
import pandas as pd
import os.path
import os
import copy
import shutil, errno
import statistics
import cometspy as c
from cobra import Reaction
from cobra import Metabolite
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import numpy as np
from pathlib import Path

############################
####FUNCIONES A UTILZAR#####
############################

#Devuelve una lista con los nombres de los modelos que crecen en un medio mínimo con los azúcares especificados
#Procede del resultado de un análisis anterior
def candidate_models(file_path):
    with open(file_path, 'r') as archivo:
        lineas = archivo.readlines()
    candidates = []
    for line in lineas:
        model_name = line.replace("\n","")
        candidates.append(model_name)
    return candidates

#La siguiente función descarga el modelo .sbml y lo guarda en una carpeta "modelos_agora", tras lo cual lo carga
#Si el modelo ya está descargado en dicha carpeta, simplemente lo carga
def descargar_leer_agora2(nombre_modelo):
    nombre_modelo = nombre_modelo
    url_padre = "https://www.vmh.life/files/reconstructions/AGORA2/version2.01/sbml_files/individual_reconstructions/"
    url = url_padre + nombre_modelo
    #Directorio donde se guardan los modelos
    nombre_directorio = "modelos_agora"
    directorio_actual = os.getcwd()
    ruta_directorio = os.path.join(directorio_actual, nombre_directorio)
    #Lo creamos si no existe
    if not os.path.exists(ruta_directorio):
        os.makedirs(ruta_directorio)
        models_dir = ruta_directorio
    else:
        models_dir = ruta_directorio
    #Si ya existe el archivo no hace falta descargarlo, se importa directamente
    ruta_archivo = os.path.join(str(models_dir), nombre_modelo)
    if os.path.exists(ruta_archivo):
        model = read_sbml_model(str(ruta_archivo))
        return model 
    #Si no existe se descarga y se importa
    else:
        try:
            # Realizar la solicitud GET a la URL del archivo
            respuesta = requests.get(url)

            # Verificar si la solicitud fue exitosa (código de estado 200)
            if respuesta.status_code == 200:
                # Guardar el archivo localmente
                with open(ruta_archivo, 'wb') as archivo_local:
                    archivo_local.write(respuesta.content)
                #print(f"Archivo descargado exitosamente como '{nombre_modelo}'.")
                #Determinamos el path al archivo
                model = read_sbml_model(str(ruta_archivo))
                return model

            else:
                print(f"Error al descargar el archivo. Código de estado: {respuesta.status_code}")

        except Exception as e:
            print("Error del modelo " + str(nombre_modelo))

################################################################
################################################################
def change_infite_value_comets(cometsfile):
    #read input file
    fin = open(cometsfile, "rt")
    #read file contents to string
    data = fin.read()
    #replace all occurrences of the required string
    data = data.replace('inf', '1000')
    #close the input file
    fin.close()
    #open the input file in write mode
    fin = open(cometsfile, "wt")
    #overrite the input file with the resulting data
    fin.write(data)
    #close the file
    fin.close()
        
def initialize_layout(models):
    #The layout is a file with the stablished format
    layout=c.layout()
    for model in models:
        layout.add_model(model)
    layout.add_typical_trace_metabolites()
    layout.update_models()
    return(layout)

def intersection_exchange_metabolites(models):
    # Calculo del medio de varias especies
    all_exchanged_mets = []
    for m in models:
        all_exchanged_mets.extend(m.get_exchange_metabolites())
    all_exchanged_mets = sorted(list(set(list(all_exchanged_mets))))
    return(all_exchanged_mets)

def create_complete_media(models):
    metabolite_list=intersection_exchange_metabolites(models)
    media_dict=dict()
    # By default 0 ammount
    for metabolite in metabolite_list:
        media_dict[metabolite]=10.0
    media_dict["ac_e"]=0.0
    media_dict["o2_e"]=1000.0
    media_dict["nh3_e"]=30.0
    media_dict["ca2_e"]=1000
    media_dict["co2_e"]=100
    media_dict["cobalt2_e"]=1000
    media_dict["cu2_e"]=1000
    media_dict["fe2_e"]=100
    media_dict["fe3_e"]=100
    media_dict["h_e"]=100
    media_dict["h2o_e"]=1000
    media_dict["k_e"]=100
    media_dict["mg2_e"]=100
    media_dict["mn2_e"]=100
    media_dict["mobd_e"]=100
    media_dict["na1_e"]=100
    media_dict["no3_e"]=100
    media_dict["photon_e"]=1000
    media_dict["pi_e"]=100
    media_dict["so4_e"]=100
    media_dict["zn2_e"]=1000
    media_dict["btn_e"]=40
    media_dict["cl_e"]=1000
    media_dict["nh4_e"]=40
    media_dict["ni2_e"]=1000

    # Return pandas dataframe
    return(pd.DataFrame(media_dict.items(),columns=['metabolite','init_amount']))

## Función que Devuelve una lista de metabolitos de las fuentes de carbono permitidas consumidas por al menos un modelo
def carbon_sources(models, sources=["glucose", "rhamnose", "xylose"]):
    available_sources = []
    for model in models:
        for rxn_id in model.medium.keys():
            rxn = model.reactions.get_by_id(rxn_id)
            met_id = list(rxn.metabolites.keys())[0].id
            met_name = list(rxn.metabolites.keys())[0].name.lower()
            for sugar in sources:
                if sugar in met_name:
                    available_sources.append(met_id)
    return(list(set(available_sources)))

#Función para traducir de "exchange" a metabolito:
def traduc_exchange(exchange):
    met = exchange.replace("EX_", "").replace("(", "[").replace(")","]")
    return met

# Funcion Para obtener un diccionario con:
# - Reaccion de exchange del medio realmente mínimo (fuentes de carbono y auxotrofías): bound
def media_cobra_dicc(model):
    #Determinamos exchange generales
    general_ex = {}
    with open('exchanges.txt', 'r') as archivo:
        lineas = archivo.readlines()

    for line in lineas:
        id, boundary = line.replace("\n", "").split(" ")
        boundary = float(boundary)
        general_ex[id] = boundary
    #Nos quedamos con las exchange del medio mínimo de cobra
    minimal_ex = list(minimal_medium(model).index)
    #Guardamos todas los ids de las reacciones del modelo
    model_rxns = [rxn.id for rxn in model.reactions]
    
    true_minimal = {}
    #Hacemos las pruebas
    with model as m:
        #Cerramos todas las exchange
        media = {}
        for id in model.medium.keys():
            media[id] = 0
        model.medium = media
        #Abrimos las del medio mínimo
        media = {}
        for reac in minimal_ex:
            media[reac] = 1000
        #Abrimos las de los metabolitos generales
        #También así se ajustan los bounds de estos exchanges
        for reac, bound in general_ex.items():
            if reac in model_rxns:
                media[reac] = bound
        #Creamos el medio
        model.medium = media

        #Las cerramos de una en una viendo si crece, si no crece lo añadimos a la lista de metabolitos minimos
        minimal_rxs = []
        for id in minimal_ex:
            with model as m:
                #Las cerramos de una en una viendo si crece
                model.reactions.get_by_id(id).bounds = (0,0)
                #print(model.medium)
                sol = model.optimize()
                if sol.objective_value < 0.1:
                    minimal_rxs.append(id)
    #Normalizamos los boundaries según número de carbono
    uptake = 5
    flujo = uptake * 6
    for rxn_id in minimal_rxs:
        try:
            n_carbonos = list(model.reactions.get_by_id(rxn_id).metabolites.keys())[0].elements["C"]
            bound = flujo/n_carbonos
            if bound > uptake:
                bound = uptake
            true_minimal[rxn_id] = bound/10 #Lo normalizamos a 10 veces menos para que no se utilice como fuente de carbono principal
        except:
            true_minimal[rxn_id] = 1000
    #Añadimos los metabolitos generales
    for reac, bound in general_ex.items():
        if reac in model_rxns:
            true_minimal[reac] = bound
 
    return true_minimal

#Ahora concatenamos todos los medios de los modelos
def general_minimal_media(models):
    general_media = {}
    for modelo in models:
        minimal_media = media_cobra_dicc(modelo)
        for exchange, bound in minimal_media.items():
            #Añadimos el exchange al medio si no está
            if exchange not in list(general_media.keys()):
                general_media[exchange] = bound
    
    return general_media

#Función para definir la concentración de azúcares en el medio dada el porcentaje respectivo a la glucosa
#Devuelve un diccionario con estos valores
#Ajustamos la concentración total a 20 mM
def sugar_levels(glc_D=18.4, xyl_D=11.6, rmn=39.4, gal=0, mnl=0, conc=20):
    sugar_names = ["glc_D", "xyl_D", "rmn", "gal", "mnl"]
    concentraciones = {}
    for name in sugar_names:
        if eval(name) > 0:
            concentraciones[name+"[e]"] = eval(name)
    #Si no suman 20 ajustamos
    suma = sum(list(concentraciones.values()))
    ratio = conc/suma
    if suma != 20:
       for azucar, conc in concentraciones.items():
           concentraciones[azucar] = round(conc*ratio, 6)
    
    return concentraciones

#Función para cerrar todos los exchange de un medio:
def close_all_ex(model):
    media_closed = {}
    for reac in model.medium.keys():
        media_closed[reac] = 0
    model.medium = media_closed

#Funcion para abrir exchanges si existen, dado un diccionario con: reaction_id, bound
def open_available_exchanges(model, dicc_exchanges):
    for reac, bound in dicc_exchanges.items():
        if reac in [rx.id for rx in model.reactions]:
            model.reactions.get_by_id(reac).lower_bound = -bound

# Función para crear un layout de COMETS con las concentraciones de matbolitos adecuadas
# Añadimos también los modelos, ajustando sus bounds
def create_complete_layout(models, biomasses, medio, ly):
    #Añadimos el medio
    medio_cobra = medio
    media_dict = {}
    for ex_id, bound in medio_cobra.items():
        met = traduc_exchange(ex_id)
        media_dict[met] = 1.5
    #Actualizamos el valor de los azúcares
    for met_id, value in sugar_levels().items():
        media_dict[met_id] = value
    
    for met, cant in media_dict.items():
        ly.set_specific_metabolite(met, cant)
    #Añadimos el peptidoglicano
    ly.set_specific_metabolite("PGPm1[c]", 100)
    #print(ly.media)
    
    #Añadimos los modelos
    for model, biomass in zip(models, biomasses):
        new_model = None
        with model as m:
            #Cerramos las exchange
            close_all_ex(model)
            #Abrimos los exchange que existan con el bound adecuado
            open_available_exchanges(model, medio_cobra)
            new_model = model.copy()
        model_added = c.model(new_model)
        model_added.initial_pop=[0,0,biomass]
        ly.add_model(model_added)

def initialize_params(package, globall):
    """Function to initialize the comets' params class
    it can be initialize by submitting two files, one for the package parameters
    and one for the global ones.
    If you don't submit a file, the params class will be initialize with the values stated below
    which have been tested in this simulation"""
    
    if package and globall is not None:
        params = c.params(global_params = globall, package_params= package)
    elif package is not None:
        params = c.params(package_params=package)
    elif globall is not None:
        params = c.params(global_params=globall)
        
    else:
        params = c.params()
        params.all_params['maxCycles']=500
        params.all_params['timeStep']=0.1
        params.all_params['spaceWidth']=0.05
        params.all_params['allowCellOverlap']= True
        params.all_params['cellSize'] = 4.3e-13
        params.all_params['deathRate']= 0.00
        params.all_params['dilFactor'] = 1e-2
        params.all_params['dilTime'] = 12
        params.all_params['numRunThreads']= 8
        params.all_params['maxSpaceBiomass']= 1000
        params.all_params['defaultVmax']= 20
        params.all_params['defaultKm']= 0.01
        params.all_params['defaultHill']= 1
        params.all_params['showCycleTime']= True
        params.all_params['geneFractionalCost']= 0
        params.all_params['writeTotalBiomassLog']= True
        params.all_params['writeMediaLog']= True
        params.all_params['writeFluxLog']= True
        params.all_params['useLogNameTimeStamp']= False
        params.all_params['flowDiffRate']= 1e-7
        params.all_params['growthDiffRate']=1e-7
        params.all_params['FluxLogRate']= 1
        params.all_params['MediaLogRate']= 1
        params.all_params['numExRxnSubsteps']= 12
        params.all_params['geneFractionalCost'] = 0
        params.all_params['exchangestyle']= 'Standard FBA'
        params.all_params['biomassMotionStyle']= 'Diffusion 2D(Crank-Nicolson)'
        params.all_params['mutRate'] = 1e-9
        params.all_params['addRate'] = 1e-9
        params.all_params['minSpaceBiomass'] = 1e-10

        
    #check if param 'maxSpaceBiomass' has default value
    if (params.all_params['maxSpaceBiomass']== 0.1):
        print('The parameter "maxSpaceBiomass" is set to 0.1.\n' \
              'It may need to change if the mo used growths well.')

    return params

def find_end_cycle(simulation_output):
    # it finish when the strains stop the growth
    end_cycle=0
    counter=0
    for index,row in simulation_output.iterrows():
        new_biomasses=row.iloc[1:]
        if index==0 or index==1:
            old_biomasses=row.iloc[1:]
        else:
            result=new_biomasses.subtract(old_biomasses)
            is_growing=result.apply(lambda x: x>0)
            old_biomasses=new_biomasses
            if is_growing.any()==False:
                counter=counter+1
                if counter==10:
                    break
        end_cycle=index
    return(end_cycle)

#Función para generar:
#1. Un .txt que contiene para cada ciclo (hasta el máximo) la biomasa de las cepas y la concentración de los metabolitos de interés
#2. Una representación gráfica de lo anterior en .pdf (cada ciclo se considera 1/10 de hora)
def make_df_and_graph(strains, metabolites, comets, max_cycles):
    '''This function creates a figure and saves it to pdf format.
    It also creates the file biomass_vs_met.txt which contais the quantity
    of each strain and metabolite and has the following columns:
    time(h), strain1 ... strainX, met1 ... metX.'''
    file_name='_'.join(metabolites)
    df = comets.media #We get the media composition results'
    df_media=copy.deepcopy(df.loc[df['cycle']<max_cycles])
    df2=comets.total_biomass #We get the biomass results
    df_biomass=copy.deepcopy(df2.loc[df2['cycle']<max_cycles])
    columns=['cycle']
    for i in range(0,len(strains)):
        columns.append(strains[i])
    df_biomass.columns=columns
    
    """For each metabolite a column with all zeros is added to the dataframe and each row that contains a value
     (metabolite concentration)is changed in the dataframe"""
    for d in metabolites:
        columns.append(d)
        met =df_media.loc[df_media['metabolite'] == d]

        temp=np.zeros(max_cycles) #Create an array with all zeros
        df_biomass[d]=temp #We added it to the dataframe
        j=1
        while j < (max_cycles): #For each cycle
            if (met.cycle==j).any(): #If the row exists
                df_biomass.loc[j-1,d] = float(met[met['cycle']==j]['conc_mmol']) #Its dataframe value is changed
            j+=1
    df_biomass.columns=columns
    np.savetxt(r'biomass_vs_'+file_name+'_template.txt', df_biomass.values, fmt='%s',delimiter='\t',header='\t'.join(columns)) #The data is saved
    
    #---Starting the figure 
    plt.ioff()
    fig, ax = plt.subplots()

    ax.set_xlabel('Tiempo (h)')
    ax.set_ylabel('Biomasa (g/L)')
    c=['r', 'g', 'b', 'orange', 'pink', 'grey']
    j=0
    lines_strains = []
    for i in strains:
        lines_strains.append(ax.plot(df_biomass['cycle']*0.1, df_biomass[i], label=" ".join(i.split("_")[:2]), color=c[j]))
        j+=1
    lines_strains = [elemento for sublist in lines_strains for elemento in sublist] #Unlist
    ax2 = ax.twinx()
    ax2.set_ylabel('Concentración de azúcares (mM)')

    #Dictionary to change sugar id to name
    sugar_id_name = {
        "rmn[e]":"Ramnosa",
        "xyl_D[e]":"Xilosa",
        "glc_D[e]":"Glucosa"
    }

    lines_met = []
    for m in metabolites:
        lines_met.append(ax2.plot(df_biomass['cycle']*0.1, df_biomass[m], label=sugar_id_name[m], color=c[j], linestyle='--'))
        j+=1
    lines_met = [elemento for sublist in lines_met for elemento in sublist] #Unlist

    #Leyendas
    font_italic = FontProperties()
    font_italic.set_style('italic')

    font_normal = FontProperties()
    font_normal.set_style('normal')

    # Crear una sola leyenda combinada
    lines = lines_strains + lines_met
    len_lines = len(lines)
    labels = [line.get_label() for line in lines]

    # Crear leyenda general
    legend = plt.legend(lines, labels, loc='upper right', prop=font_normal, ncol=1)

    # Obtener las etiquetas de la leyenda y asignarles el estilo de fuente adecuado
    for text, font in zip(legend.get_texts(), [font_italic] * 3 + [font_normal] * (len_lines-3)):
        text.set_fontproperties(font)
    
    #saving the figure to a pdf
    plt.savefig('biomass_vs_'+file_name+'_template_plot.pdf')
    
    return df_biomass

#####################
########FLYCOP#######
#####################

#Set matplotlib backend for batch operation
import matplotlib
matplotlib.use("Agg")

#Lista de organismos a iterar
lista_organismos = candidate_models("selected_models.txt")

#Medio general definido anteriormente
import json
# Abrir el archivo JSON y cargar el contenido en un diccionario
with open('Medio_RQ.json', 'r') as archivo:
    medio_general = json.load(archivo)

def FLYCOP_space():
    # Build Configuration Space which defines all parameters and their ranges
    cs = ConfigurationSpace()
    e1 = Categorical('organism_1', lista_organismos, default="Bacillus_amyloliquefaciens_subsp_amyloliquefaciens_DC_12.xml")
    e2 = Categorical('organism_2', lista_organismos, default="Bacillus_amyloliquefaciens_subsp_amyloliquefaciens_DC_12.xml")
    e3 = Categorical('organism_3', lista_organismos, default="Bacillus_amyloliquefaciens_subsp_amyloliquefaciens_DC_12.xml")
    b1 = Float("biomass_1", bounds=(0.0001,0.1), default=0.01)
    b2 = Float("biomass_2", bounds =(0.0001,0.1), default=0.01)
    b3 = Float("biomass_3", bounds =(0.0001,0.1), default=0.01)
    o1 = Constant("output_file",value="promiseang_bottom_up.csv")
    
    cs.add_hyperparameters([e1,e2,e3,b1,b2,b3,o1])

    """
    #Conditions for not repeating organisms
    for nombre in lista_organismos:
        cs.add_forbidden_clause(ForbiddenAndConjunction(
            ForbiddenEqualsClause(e1, nombre),
            ForbiddenEqualsClause(e2, nombre)
        ))
        cs.add_forbidden_clause(ForbiddenAndConjunction(
            ForbiddenEqualsClause(e1, nombre),
            ForbiddenEqualsClause(e3, nombre)
        ))
        # Prohibir que nombre2 sea igual a nombre3
        cs.add_forbidden_clause(ForbiddenAndConjunction(
            ForbiddenEqualsClause(e2, nombre),
            ForbiddenEqualsClause(e3, nombre)
        ))
    """
    
    scenario = Scenario(
        cs,
        deterministic = True,
        output_directory=Path("./PROMISEANG_FLYCOP_results"),
        #walltime_limit=120,  # Limit to two minutes
        n_trials=1000,  # Evaluated max 500 trials
        #n_workers=8,  # Use eight workers
    )

    return(scenario)

def print_scenario(parameters):
    for key, value in parameters.items():
        print(key, "->", value)

def FLYCOP_schema(config: Configuration, seed: int = 0) -> float:
    #load parameters
    e1 = config["organism_1"]
    e2 = config["organism_2"]
    e3 = config["organism_3"]
    b1 = config["biomass_1"]
    b2 = config["biomass_2"]
    b3 = config["biomass_3"]
    o1 = config["output_file"]
 
    parameters={}
    parameters["organism_1"] = e1
    parameters["organism_2"] = e2
    parameters["organism_3"] = e3
    parameters["biomass_1"] = b1
    parameters["biomass_2"] = b2
    parameters["biomass_3"] = b3
    parameters["output_file"] = o1

    print_scenario(parameters)
    
    #Biomasses
    biomass1=parameters["biomass_1"]
    biomass2=parameters["biomass_2"]
    biomass3=parameters["biomass_3"]
    
    #Strains
    strain1 = parameters["organism_1"].replace(".xml","")
    strain2 = parameters["organism_2"].replace(".xml","")
    strain3 = parameters["organism_3"].replace(".xml","")
    strains=[strain1, strain2, strain3]
    
    #####################
    # COMETS run        #
    #####################
    
    # COBRA MODELS
    bacteria_1_cobra = descargar_leer_agora2(parameters["organism_1"])
    bacteria_2_cobra = descargar_leer_agora2(parameters["organism_2"])
    bacteria_3_cobra = descargar_leer_agora2(parameters["organism_3"])
 
    # Layout
    layout = c.layout()
    create_complete_layout([bacteria_1_cobra, bacteria_2_cobra, bacteria_3_cobra],[biomass1, biomass2, biomass3], medio_general, layout)

    # Metabolites to analyse (available sugars which can be consumed by the bacteria)
    metabolites = carbon_sources([bacteria_1_cobra, bacteria_2_cobra, bacteria_3_cobra])

    print("Loading the comets parameters...")
    params_package=None
    params_global=None
    params =initialize_params(params_package, params_global)
    print("Creating the comets object...")
    comets = c.comets(layout, params)
    dirPlot='promiseang_output'
    comets.run(delete_files=True)
    
    comets.total_biomass.to_csv('Total_biomass_log_template.txt',sep='\t',index=False)
    endCycle=find_end_cycle(comets.total_biomass)
        
    print("The last cycle is: %d"%endCycle)
    maxCycles=endCycle-1
    
    # Plot
    df2=make_df_and_graph(strains, metabolites, comets, endCycle)
    final_biomass1=float(df2.at[maxCycles-1,strain1])
    final_biomass2=float(df2.at[maxCycles-1,strain2])
    final_biomass3=float(df2.at[maxCycles-1,strain3])
    final_biomasses = [final_biomass1, final_biomass2, final_biomass3]

    # Fitness(penalizes big differences among biomasses)
    # The consortia must consume all sugars of interest
    total_biomass = sum(final_biomasses)
    dispersion = max(final_biomasses)-min(final_biomasses)
    if len(metabolites) < 3:
        fitness = 0
    else:
        fitness = total_biomass - dispersion

    # If two species are the same, fitness is 0
    if len({strain1, strain2, strain3})<3:
        fitness = 0

    #If after the penalties application, fitness is negative, change it to zero:
    if fitness < 0:
        fitness = 0

    # Create a string with the parameters separated by comma
    configuration=','.join([str(value)for value in parameters.values()])
    sugars_consumed='\t'.join([str(metabolite)for metabolite in metabolites])
    header_line="Configuration\tfitness\t"+"Sugars_consumed\t"+"\tBiomass1\tBiomass2\tBiomass3\tEndCycle"
    print(header_line)
    output_line=configuration+'\t'+str(fitness)+'\t'+sugars_consumed+'\t'+str(final_biomass1)+'\t'+str(final_biomass2)+'\t'+str(final_biomass3)+'\t'+str(endCycle)
    print(output_line)
    if exists(parameters["output_file"]):
        with open(parameters["output_file"], 'a') as file:
            file.write(output_line+"\n")
    else:
        with open(parameters["output_file"], 'a') as file:
            file.write(header_line+"\n")
            file.write(output_line+"\n")
            
    fitness_inverse=-1*fitness
    return fitness_inverse,0

def FLYCOP_run():
    print("Optimizing! Depending on your machine, this might take a few minutes.")
    scenario = FLYCOP_space()
    #print(scenario)
    smac = HyperparameterOptimizationFacade(scenario, FLYCOP_schema)
    result=smac.optimize()
    return(result)

##################
#####PRUEBA#######
##################
"""
parameters={
    "organism_1":'Corynebacterium_casei_UCMA_3821.xml',
    "organism_2":'Bacillus_paralicheniformis_ATCC_9945a.xml',
    "organism_3":'Bacillus_amyloliquefaciens_subsp_amyloliquefaciens_DC_12.xml',
    "biomass_1":0.04,
    "biomass_2":0.03,
    "biomass_3":0.01,
    "output_file":'promiseang_bottom_up.csv'}
FLYCOP_schema(parameters)
"""

FLYCOP_run()


