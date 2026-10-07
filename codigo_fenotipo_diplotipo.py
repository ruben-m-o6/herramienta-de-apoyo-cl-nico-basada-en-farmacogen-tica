import pandas as pd

# 1. Cargar el archivo CSV original
df = pd.read_csv("Genotype Matrix.csv", sep=";")
df.columns = [col.strip() for col in df.columns]

# --- 2. FUNCIONES PARA EVALUAR CADA VARIANTE INDIVIDUALMENTE ---

# GEN CYP2D6
def evaluar_cyp2d6_3(val):
    if val == 'T/T': return '*1/*1'      
    elif val == 'T/-': return '*1/*3'    
    else: return '*3/*3'                 

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
    if val == 'C/C': return '*1/*1'
    elif val == 'C/T': return '*1/*28'
    elif val == 'T/T': return '*28/*28'
    else: return '*1/*1'

# --- 3. APLICAR FUNCIONES INDIVIDUALES ---
df['Genotipo_CYP2D6*3'] = df['CYP2D6*3'].apply(evaluar_cyp2d6_3)
df['Genotipo_CYP2D6*4'] = df['CYP2D6*4'].apply(evaluar_cyp2d6_4)
df['Genotipo_CYP2D6*17'] = df['CYP2D6*17'].apply(evaluar_cyp2d6_17)
df['Genotipo_CYP2D6*41'] = df['CYP2D6*41'].apply(evaluar_cyp2d6_41)

df['Genotipo_DPYD*2A'] = df['DPYD*2A'].apply(evaluar_dpyd_2A)
df['Genotipo_DPYD*13'] = df['DPYD*13'].apply(evaluar_dpyd_13)
df['Genotipo_DPYD_HapB3'] = df['DPYD_HapB3'].apply(evaluar_dpyd_HapB3)
df['Genotipo_DPYD_D949V'] = df['DPYD_D949V'].apply(evaluar_dpyd_D949V)

df['Genotipo_UGT1A1*80'] = df['UGT1A1*80'].apply(evaluar_ugt1a1_80)


# --- 4. NUEVO PASO: CALCULAR DIPLOTIPO (Integrar mutaciones) ---
def obtener_diplotipo(row, columnas_gen):
    """Extrae todos los alelos mutados de un gen y los fusiona en un solo diplotipo"""
    alelos_mutados = []
    for col in columnas_gen:
        alelos = row[col].split('/')
        for alelo in alelos:
            if alelo != '*1': # Ignoramos el *1 porque es el "normal"
                alelos_mutados.append(alelo)
    
    # Regla de Oro: Si no hay mutaciones, es *1/*1
    if len(alelos_mutados) == 0:
        return '*1/*1'
    elif len(alelos_mutados) == 1:
        return f'*1/{alelos_mutados[0]}'
    elif len(alelos_mutados) == 2:
        return f'{alelos_mutados[0]}/{alelos_mutados[1]}'
    else:
        # Por si hay más de 2 mutaciones (casos muy raros)
        return f'{alelos_mutados[0]}/{alelos_mutados[1]}'

cols_cyp2d6 = ['Genotipo_CYP2D6*3', 'Genotipo_CYP2D6*4', 'Genotipo_CYP2D6*17', 'Genotipo_CYP2D6*41']
cols_dpyd = ['Genotipo_DPYD*2A', 'Genotipo_DPYD*13', 'Genotipo_DPYD_HapB3', 'Genotipo_DPYD_D949V']

df['Diplotipo_CYP2D6'] = df.apply(lambda row: obtener_diplotipo(row, cols_cyp2d6), axis=1)
df['Diplotipo_DPYD'] = df.apply(lambda row: obtener_diplotipo(row, cols_dpyd), axis=1)
df['Diplotipo_UGT1A1'] = df['Genotipo_UGT1A1*80'] # Como solo evaluamos un SNP, el genotipo individual es directamente el diplotipo


# --- 5. NUEVO PASO: CALCULAR FENOTIPO (Score de Actividad) ---

def fenotipo_ugt1a1(diplotipo):
    if diplotipo == '*1/*1': return 'NM (Metabolizador Normal)'
    elif diplotipo == '*1/*28': return 'IM (Metabolizador Intermedio)'
    elif diplotipo == '*28/*28': return 'PM (Metabolizador Lento)'
    return 'Indeterminado'

def fenotipo_dpyd(diplotipo):
    # Valores según CPIC: Normal=1, Deficientes=0, Reducidas=0.5
    scores = {'*1': 1.0, '*2A': 0.0, '*13': 0.0, 'HapB3': 0.5, 'D949V': 0.5}
    alelos = diplotipo.split('/')
    # Sumar el score de los 2 alelos
    score_total = sum(scores.get(a, 1.0) for a in alelos)
    
    if score_total == 2.0: return 'NM (Metabolizador Normal)'
    elif 1.0 <= score_total <= 1.5: return 'IM (Metabolizador Intermedio)'
    elif score_total <= 0.5: return 'PM (Metabolizador Lento)'
    return 'Indeterminado'

def fenotipo_cyp2d6(diplotipo):
    # Valores según CPIC: Normal=1, Nulos=0, Reducidas=0.5
    scores = {'*1': 1.0, '*3': 0.0, '*4': 0.0, '*17': 0.5, '*41': 0.5}
    alelos = diplotipo.split('/')
    score_total = sum(scores.get(a, 1.0) for a in alelos)
    
    if 1.25 <= score_total <= 2.25: return 'NM (Metabolizador Normal)'
    elif 0.25 <= score_total <= 1.0: return 'IM (Metabolizador Intermedio)'
    elif score_total == 0.0: return 'PM (Metabolizador Lento)'
    # Nota de tus apuntes: con TaqMan genérico no solemos detectar UM (duplicaciones), pero lo incluimos
    elif score_total > 2.25: return 'UM (Metabolizador Ultrarrápido)' 
    return 'Indeterminado'

df['Fenotipo_CYP2D6'] = df['Diplotipo_CYP2D6'].apply(fenotipo_cyp2d6)
df['Fenotipo_DPYD'] = df['Diplotipo_DPYD'].apply(fenotipo_dpyd)
df['Fenotipo_UGT1A1'] = df['Diplotipo_UGT1A1'].apply(fenotipo_ugt1a1)

# 6. Guardar la tabla extendida con los resultados listos para generar el informe
columnas_finales = ['Sample/Assay', 'Diplotipo_CYP2D6', 'Fenotipo_CYP2D6', 
                    'Diplotipo_DPYD', 'Fenotipo_DPYD', 
                    'Diplotipo_UGT1A1', 'Fenotipo_UGT1A1']

df_final = df[columnas_finales]
df_final.to_csv("Matriz_Resultados_Completos.csv", index=False, sep=";")

print("¡Proceso exitoso! Se han calculado los Diplotipos y Fenotipos.")
print(df_final.head())