import pandas as pd
from fpdf import FPDF
import os
import tkinter as tk
from tkinter import filedialog, messagebox

def procesar_archivo():
    filepath = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
    if not filepath:
        return
    
    try:
        # 1. Leer y limpiar el CSV original
        df = pd.read_csv(filepath, sep=";")
        df.columns = [col.strip() for col in df.columns]
        
        # (Aquí va la lógica de los diplotipos y fenotipos que ya programaste)
        # Para mantener este código corto, asumimos que aquí cruzas los datos
        # y generas el dataframe final "df_final" listo para los PDFs.
        
        # 2. Generar carpeta de salida
        if not os.path.exists("Informes_Clinicos"):
            os.makedirs("Informes_Clinicos")
            
        # 3. Simulación rápida de la generación para la interfaz
        # (Asegúrate de pegar aquí las funciones generar_pdf() del paso anterior)
        
        messagebox.showinfo("Éxito", f"Se han procesado {len(df)} pacientes y generado los PDFs en la carpeta 'Informes_Clinicos'.")
        
    except Exception as e:
        messagebox.showerror("Error", f"Ha ocurrido un error: {str(e)}")

# --- INTERFAZ GRÁFICA ---
root = tk.Tk()
root.title("Herramienta Farmacogenética - Oncología")
root.geometry("400x200")

label = tk.Label(root, text="Generador de Informes Clínicos", font=("Arial", 14, "bold"))
label.pack(pady=20)

btn_cargar = tk.Button(root, text="Cargar CSV y Generar PDFs", command=procesar_archivo, bg="lightblue", font=("Arial", 12))
btn_cargar.pack(pady=10)

root.mainloop()
