import pandas as pd

# 1. Cargar el archivo CSV original
df = pd.read_csv("Genotype Matrix.csv", sep=";")
df.columns = [col.strip() for col in df.columns]

# --- 2. FUNCIONES PARA EVALUAR CADA VARIANTE INDIVIDUALMENTE ---

# GEN CYP2D6
def evaluar_cyp2d6_3(val):
    if val == 'T/T': return '*1/*1'      # Ancestral (Normal)
    elif val == 'T/-': return '*1/*3'    # Heterocigoto
    else: return '*3/*3'                 # Homocigoto mutado

def evaluar_cyp2d6_4(val):
    if val == 'C/C': return '*1/*1'
    elif val == 'C/T': return '*1/*4'
    elif val == 'T/T': return '*4/*4'
    else: return '*1/*1'

def evaluar_cyp2d6_17(val):
    if val == 'G/G': return '*1/*1'
    elif val == 'G/A': return '*1/*17'
    elif val == 'A/A': return '*17/*17'
    else: return '*1/*1'

def evaluar_cyp2d6_41(val):
    if val == 'C/C': return '*1/*1'
    elif val == 'C/T': return '*1/*41'
    elif val == 'T/T': return '*41/*41'
    else: return '*1/*1'

# GEN DPYD
def evaluar_dpyd_2A(val):
    if val == 'C/C': return '*1/*1'
    elif val == 'C/T' or val == 'G/C': return '*1/*2A' 
    else: return '*2A/*2A'

def evaluar_dpyd_13(val):
    if val == 'A/A': return '*1/*1'
    elif val != 'A/A' and val[0] != val[2]: return '*1/*13'
    else: return '*13/*13'

def evaluar_dpyd_HapB3(val):
    if val == 'C/C': return '*1/*1'
    elif val == 'C/T': return '*1/HapB3'
    else: return 'HapB3/HapB3'

def evaluar_dpyd_D949V(val):
    if val == 'T/T': return '*1/*1'
    elif val != 'T/T' and val[0] != val[2]: return '*1/D949V'
    else: return 'D949V/D949V'

# GEN UGT1A1
def evaluar_ugt1a1_80(val):
    # *80 es marcador del alelo *28
    if val == 'C/C': return '*1/*1'
    elif val == 'C/T': return '*1/*28'
    elif val == 'T/T': return '*28/*28'
    else: return '*1/*1'


# --- 3. APLICAR FUNCIONES PARA CREAR LAS NUEVAS COLUMNAS ---

# Creamos una nueva columna de "Genotipo" por cada variante leída
df['Genotipo_CYP2D6*3'] = df['CYP2D6*3'].apply(evaluar_cyp2d6_3)
df['Genotipo_CYP2D6*4'] = df['CYP2D6*4'].apply(evaluar_cyp2d6_4)
df['Genotipo_CYP2D6*17'] = df['CYP2D6*17'].apply(evaluar_cyp2d6_17)
df['Genotipo_CYP2D6*41'] = df['CYP2D6*41'].apply(evaluar_cyp2d6_41)

df['Genotipo_DPYD*2A'] = df['DPYD*2A'].apply(evaluar_dpyd_2A)
df['Genotipo_DPYD*13'] = df['DPYD*13'].apply(evaluar_dpyd_13)
df['Genotipo_DPYD_HapB3'] = df['DPYD_HapB3'].apply(evaluar_dpyd_HapB3)
df['Genotipo_DPYD_D949V'] = df['DPYD_D949V'].apply(evaluar_dpyd_D949V)

df['Genotipo_UGT1A1*80'] = df['UGT1A1*80'].apply(evaluar_ugt1a1_80)

# 4. Guardar la tabla extendida
df.to_csv("Matriz_con_Genotipos_Individuales.csv", index=False, sep=";")

print("Tabla generada correctamente. Columnas añadidas:")
print(df.columns.tolist())