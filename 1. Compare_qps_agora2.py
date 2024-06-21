# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 11:07:20 2024

@author: rqr20
"""

import openpyxl
import requests
from bs4 import BeautifulSoup

#Extraemos los nombres de los organismos y los introducimos en una lista

def leer_excel(file_path):
    try:
        # Cargar el libro de trabajo de Excel
        workbook = openpyxl.load_workbook(file_path)
        
        # Seleccionar la hoja de trabajo (puedes cambiar el nombre según tu archivo)
        sheet = workbook.active
        
        # Lista para almacenar los datos de la primera columna
        datos_primera_columna = []

        # Iterar sobre las filas a partir de la tercera fila
        for row in sheet.iter_rows(min_row=3, max_col=1, values_only=True):
            # Agregar el contenido de la primera columna a la lista
            datos_primera_columna.append(str(row[0]))

        return datos_primera_columna

    except Exception as e:
        print(f"Error al leer el archivo Excel: {e}")
        return None

# Ruta de tu archivo Excel
archivo_excel = 'QPS_org.xlsx'

# Llamar a la función y obtener la lista de datos
lista_datos = leer_excel(archivo_excel)

#Formateamos los nombres
lista_nombres = []

for dato in lista_datos:
    if dato.endswith(" "):
        dato = dato[:-1]
    dato = dato.replace(" ", "_")
    lista_nombres.append(dato)

#Imprimimos la lista de organismos
if lista_nombres:
    print("Lista de organismos aceptados recuperada\n")
#    print('La lista de organismos aceptados es:')
#    for nombre in lista_nombres:
#        print(nombre)
else:
    print("No se pudo leer el archivo Excel.")
    

#Ahora extraemos los nombres de los modelos disponibles en AGORA2 y los formateamos
def obtener_enlaces_archivos(url):
    try:
        # Realizar la solicitud HTTP a la URL
        respuesta = requests.get(url)
        respuesta.raise_for_status()  # Verificar si la solicitud fue exitosa

        # Obtener el contenido HTML de la respuesta
        contenido_html = respuesta.text

        # Crear un objeto BeautifulSoup para analizar el contenido HTML
        soup = BeautifulSoup(contenido_html, 'html.parser')

        # Obtener todos los elementos <a> que contienen enlaces
        enlaces = soup.find_all('a', href=True)

        # Filtrar los enlaces que apuntan a archivos
        enlaces_archivos = [enlace['href'] for enlace in enlaces if es_enlace_archivo(enlace['href'])]

        return enlaces_archivos

    except requests.exceptions.RequestException as e:
        print(f"Error al obtener la página: {e}")
        return None

def es_enlace_archivo(enlace):
    # Puedes ajustar esta función según los tipos de archivos que estás buscando
    tipos_de_archivo_permitidos = ['.xml']
    return any(enlace.endswith(tipo) for tipo in tipos_de_archivo_permitidos)

# URL de la página con enlaces a archivos
url_pagina_con_archivos = 'https://www.vmh.life/files/reconstructions/AGORA2/version2.01/sbml_files/individual_reconstructions/'

# Obtener los enlaces a archivos desde la página
enlaces_archivos = obtener_enlaces_archivos(url_pagina_con_archivos)

#Obtenemos la lista de modelos disponibles
if enlaces_archivos:
    print("Lista de archivos recuperada")
#    print("Enlaces a archivos encontrados:")
#    for enlace in enlaces_archivos:
#        print(enlace)
else:
    print("No se encontraron enlaces a archivos en la página.")
    
#Creamos un diccionario con keys cada organismo aceptado y values los modelos para ese organismo
organismo_modelos = {}
organismo_soportado = {} #Solo nos dice si hay modelos o no
genero_modelos = {} # Lo mismo pero para géneros
genero_soportado = {}


#Recuperamos para cada especie la lista de modelos de cepas disponibles en AGORA2, y si está soportado o no
#En caso de que no haya modelos para la especie, buscamos para el género
for nombre in lista_nombres:
    #print(nombre)
    cepas = [] #lista de modelos para cada organismo
    for modelo in enlaces_archivos:
        #print(modelo)
        if nombre in modelo:
            cepas.append(modelo)
    #Metemos entradas a los diccionario
    organismo_modelos[nombre] = cepas
    #Comprobamos si la lista está vacía
    if not cepas:
        organismo_soportado[nombre] = "N"

    else:
        organismo_soportado[nombre] = "Y"
        
#Repetimos el procedimiento con los géneros
for nombre in lista_nombres:
    genero = nombre.split("_")[0] #Nos quedamos con el género
    modelos_gen = []
    for modelo in enlaces_archivos:
        if genero in modelo:
            modelos_gen.append(modelo)
    #Metemos entrada al diccionario        
    genero_modelos[genero] = modelos_gen
    if not modelos_gen:
        genero_soportado[genero] = "N"
        
    else:
        genero_soportado[genero] = "Y"

# Pasamos los resultados a txts
with open("QPS_especies_soportadas.txt", 'w') as archivo:
    # Escribimos
    for nombre, soportado in organismo_soportado.items() :
        archivo.write(nombre + ": " + soportado + "\n")

with open("QPS_especies_modelos.txt", 'w') as archivo:
    # Escribimos
    for nombre, modelos in organismo_modelos.items() :
        archivo.write(nombre + ": ")
        for mod in modelos:
            archivo.write(mod + ", ")
        archivo.write("\n")
        
with open("QPS_generos_soportados.txt", 'w') as archivo:
    # Escribimos
    for nombre, soportado in genero_soportado.items() :
        archivo.write(nombre + ": " + soportado + "\n")

with open("QPS_generos_modelos.txt", 'w') as archivo:
    # Escribimos
    for nombre, modelos in genero_modelos.items() :
        archivo.write(nombre + ": ")
        for mod in modelos:
            archivo.write(mod + ", ")
        archivo.write("\n")
        
#Transformamos las listas de organismos y generos soportados a formato excel
import pandas as pd

# Crear un DataFrame a partir del diccionario y añadir una columna con los modelos
df_esp = pd.DataFrame(list(organismo_soportado.items()), columns=['Especie', 'Modelo'])
df_gen = pd.DataFrame(list(genero_soportado.items()), columns=['Genero', 'Modelo'])

#Para cada organismo unimos todos los modelos que hay de él en un string, separados por "," y creamos una lista
lista_modelos_org = [','.join(sublista) for sublista in organismo_modelos.values()]
lista_modelos_gen = [','.join(sublista) for sublista in genero_modelos.values()]

#Modificamos los dataframes para añadir la tercera columna con los modelos
df_esp["Modelos"] = lista_modelos_org
df_gen["Modelos"] = lista_modelos_gen

# Especifica el nombre del archivo de Excel
excel_especies = 'especies_soportadas.xlsx'
excel_generos = 'generos_soportados.xlsx'

# Guarda el DataFrame en un archivo de Excel
df_esp.to_excel(excel_especies, index=False)
df_gen.to_excel(excel_generos, index = False)

