# Importing module
import tkinter as tk
from tkinter import filedialog, StringVar
from tkinter import messagebox
from pathlib import Path
import pandas as pd
import hashlib
from threading import Thread
import numpy as np
from PIL import Image, ImageTk, ImageOps

from gwas_analysis.plots import plot_joint_manhattan, plot_scatter_with_regression
from gwas_analysis.stats import tidy_summary_stats
from gwas_analysis.ldsc import ldsc, LDSCResult

# Defining Function
# Function to upload a file and display its content in the GUI
def upload_file(upload_label):
    file_path = filedialog.askopenfilename(filetypes=[("TSV files", "*.tsv")])
    if file_path:
        confirm = messagebox.askyesno("Confirm Upload", f"Are you sure you want to upload this file?\n{file_path}")
        if confirm:
            if upload_label.winfo_name() == "disease1":
                disease1_path.set(file_path)
            else:
                disease2_path.set(file_path)
            upload_label.config(text=f"{shorten_path(file_path)}", font = ("Helvetica", round(user_height/70)))
            check_files_uploaded()
        else:
            upload_label.config(text="Upload Cancelled")
            upload_label.config(text="No File Selected")

def shorten_path(full_path):
    p = Path(full_path)
    q = p.parts[len(p.parts)-1]
    if len(p.parts) > 2:
        if len(q) >= 20:
            r = q[0:6]+"..."+q[len(q)-10:len(q)]
            return ".../"+r
        else:
            return ".../"+q
    else:
        return str(p)

def load_data():
    if stats_test.get() == "LDSC":
        root.after(0, lambda: go.config(text="Calculating LDSC..."))
        result = ldsc(disease1_path.get(), ldsc_n1.get(), disease2_path.get(), ldsc_n2.get(), ldsc_stat1.get(), ldsc_stat2.get())
        root.after(0, display_ldsc, result)
        return

    with open(disease1_path.get(), 'rb', buffering=0) as f:
        hash1 = hashlib.file_digest(f, 'md5').hexdigest()
    with open(disease2_path.get(), 'rb', buffering=0) as f:
        hash2 = hashlib.file_digest(f, 'md5').hexdigest()

    try:
        root.after(0, lambda: go.config(text="Reading merged data file"))
        merged_df = pd.read_csv(f'data/{hash1}_{hash2}.tsv', sep='\t')
    except FileNotFoundError:
        root.after(0, lambda: go.config(text="Loading disease 1"))
        disease1 = tidy_summary_stats(pd.read_csv(disease1_path.get(), sep='\t'), significance=1)
        root.after(0, lambda: go.config(text="Loading disease 2"))
        disease2 = tidy_summary_stats(pd.read_csv(disease2_path.get(), sep='\t'), significance=1)
        root.after(0, lambda: go.config(text="Merging datasets"))
        merged_df = pd.merge(disease1, disease2, how='inner', on=['risk_allele', 'chromosome', 'base_pair_location', 'effect_allele', 'other_allele'], suffixes=("_disease1", "_disease2"))
        root.after(0, lambda: go.config(text="Merging datasets"))
        merged_df.to_csv(f'data/{hash1}_{hash2}.tsv', index=False, sep='\t')

    root.after(0, lambda: go.config(text="Computing statistical method"))
    if stats_test.get() == "figures":
        filtered_df = merged_df[(merged_df['p_value_disease1'] < 1e-5) | (merged_df['p_value_disease2'] < 1e-5)]
        # filtered_df = filtered_df[(np.abs(filtered_df['beta_disease1']) > 0.01) | (np.abs(filtered_df['beta_disease2']) > 0.01)]
        plot_scatter_with_regression(filtered_df, col1='beta_disease1', col2='beta_disease2', output="beta_regression")
        plot_scatter_with_regression(filtered_df, col1='neg_log_p_value_disease1', col2='neg_log_p_value_disease2', output="p_regression")
        plot_joint_manhattan(merged_df, 'neg_log_p_value_disease1', 'neg_log_p_value_disease2')
        root.after(0, display_images)
    elif stats_test.get() == "sig":
        filtered_df = merged_df[(merged_df['p_value_disease1'] < 1e-5) | (merged_df['p_value_disease2'] < 1e-5)]
        # filtered_df: pd.DataFrame = filtered_df[(np.abs(filtered_df['beta_disease1']) > 0.05) | (np.abs(filtered_df['beta_disease2']) > 0.05)]
        filtered_df.to_csv('data/significant_SNPs.tsv', index=False, sep='\t')
        root.after(0, display_significant, filtered_df)


def display_significant(df: pd.DataFrame):
    result = """
    Significant SNPs saved to data/significant_SNPs.tsv
    """
    result_label1.config(text=result, image='')
    result_label1.grid()
    result_label2.grid_remove()
    result_label3.grid_remove()
    go.config(text="GO!", state="normal")

def display_images():
    result_image1 = tk.PhotoImage(file='data/joint_manhattan.png')
    result_label1.config(text='', image=result_image1)
    result_label1.image = result_image1
    result_label1.grid()
    size = (round(user_width*0.7*0.5), round(user_height*0.7*0.5))
    result_image2 = ImageOps.contain(Image.open('data/beta_regression.png'), size)
    resized2 = ImageTk.PhotoImage(result_image2)
    result_label2.config(text='', image=resized2)
    result_label2.image = resized2
    result_label2.grid()
    result_image3 = ImageOps.contain(Image.open('data/p_regression.png'), size)
    resized3 = ImageTk.PhotoImage(result_image3)
    result_label3.config(text='', image=resized3)
    result_label3.image = resized3
    result_label3.grid()
    go.config(text="GO!", state="normal")


def display_ldsc(result: LDSCResult | None):
    if result:
        str_result = f"""
        Genetic correlation: {result.correlation}
        Standard error: {result.correlation_std_error}
        P value: {result.p_value}
        """
    else:
        str_result = "Unable to compute LDSC for provided datasets"
    result_label1.config(text=str_result, image='')
    result_label1.grid()
    result_label2.grid_remove()
    result_label3.grid_remove()
    go.config(text="GO!", state="normal")


def GO():
    if not check_files_uploaded(): return
    go.config(text="Loading...", state="disabled")
    thread = Thread(target=load_data)
    thread.start()


def check_files_uploaded():
    if lbl_inputs_1.cget("text") != "No file selected" and lbl_inputs_2.cget("text") != "No file selected":
        go.config(bg="green")
        return True
    return False

def integer_validate(x: str):
    return x.isdigit() or not x

def show_test_variables():
    if stats_test.get() == "LDSC":
        lbl_inputs_4.grid()
        lbl_ldsc_1.grid()
        lbl_ldsc_2.grid()
        lbl_ldsc_3.grid()
        lbl_ldsc_4.grid()
        ldsc_input_1.grid()
        ldsc_input_2.grid()
        ldsc_input_3.grid()
        ldsc_input_4.grid()
    else:
        lbl_inputs_4.grid_remove()
        lbl_ldsc_1.grid_remove()
        lbl_ldsc_2.grid_remove()
        lbl_ldsc_3.grid_remove()
        lbl_ldsc_4.grid_remove()
        ldsc_input_1.grid_remove()
        ldsc_input_2.grid_remove()
        ldsc_input_3.grid_remove()
        ldsc_input_4.grid_remove()

if __name__ == "__main__":
    Path("data").mkdir(parents=True, exist_ok=True)
    # Setting up the window
    # Defining an instance window called "root" into Tk class
    root = tk.Tk()
    root.title("CLAS: Cross-Linkage Association Studies")
    root.geometry("700x700")
    root.configure(bg="#E3F2FD")
    # Getting user's window information
    user_width = root.winfo_screenwidth()
    user_height = root.winfo_screenheight()


    # Defining Widgets
    # Create frame for title: frame_title
    frame_title = tk.Frame(root, width=user_width, height=user_height/12, bg="#003366")
    # Prevent the frame from shrinking to fit its content
    frame_title.pack_propagate(False)
    # Create frame for content: frame_content
    frame_content = tk.Frame(root, width=user_width, height=(user_height-(user_height/12)))
    frame_content.grid_propagate(False)
    #Create frame for inputs: frame_inputs
    frame_inputs = tk.Frame(frame_content, width=(user_width*0.30), height=(user_height-(user_height/12)))
    frame_inputs.grid_propagate(False)
    #Create frame for outputs: frame_outputs
    frame_outputs = tk.Frame(frame_content, width=(user_width*0.70), height=(user_height-(user_height/12)))
    frame_outputs.grid_propagate(False)


    # Define title in the frame_title
    lbl_title = tk.Label(frame_title, text="CLAS: Cross-Linkage Association Studies",
                         font=("Roboto", round(user_height/25), "bold"), fg="white", bg="#003366")
    # Ensure the frame fills the available space
    frame_title.pack(fill="both", expand=True)
    # Center the label within the frame
    lbl_title.pack(expand=True)


    # Defining label in frame_inputs
    lbl_inputs_0 = tk.Label(frame_inputs, text="Upload GWAS Summary (.tsv):", font=("Roboto", round(user_height/45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)
    lbl_inputs_1 = tk.Label(frame_inputs, text="No file selected", name='disease1', font=("Roboto", round(user_height/70)))
    lbl_inputs_2 = tk.Label(frame_inputs, text="No file selected", name='disease2', font=("Roboto", round(user_height/70)))
    lbl_inputs_3 = tk.Label(frame_inputs, text="Select Method:", font=("Roboto", round(user_height/45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)

    lbl_outputs_1 = tk.Label(frame_outputs, text="Results:", font=("Roboto", round(user_height / 45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)

    # Defining button in frame_inputs
    upload_1 = tk.Button(frame_inputs, text="Genetic Disease 1", font=("Roboto", round(user_height/70)), command=lambda: upload_file(lbl_inputs_1))
    upload_2 = tk.Button(frame_inputs, text="Genetic Disease 2", font=("Roboto", round(user_height/70)), command=lambda: upload_file(lbl_inputs_2))
    go = tk.Button(frame_inputs, text="GO!", font=("Roboto", round(user_height/60), "bold"), bg="#FF3B3B", fg="black", padx=20, pady=10, command=lambda: GO())

    # Defining radio buttons in frame_inputs
    stats_test = tk.StringVar(value="figures")
    radio1 = tk.Radiobutton(frame_inputs, text="Generate Figures", font=("Roboto", round(user_height/70)), variable=stats_test, value="figures", command=show_test_variables)
    radio2 = tk.Radiobutton(frame_inputs, text="Linkage Disequilibrium Score Regression (LDSC)", font=("Roboto", round(user_height/70)), variable=stats_test, value="LDSC", command=show_test_variables)
    radio4 = tk.Radiobutton(frame_inputs, text="Significant SNPs", font=("Roboto", round(user_height/70)), variable=stats_test, value="sig", command=show_test_variables)

    lbl_inputs_4 = tk.Label(frame_inputs, text="Select Options:", font=("Roboto", round(user_height/45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)
    ldsc_n1 = tk.IntVar(value=0)
    ldsc_n2 = tk.IntVar(value=0)
    ldsc_stat1 = tk.StringVar(value='beta')
    ldsc_stat2 = tk.StringVar(value='beta')
    lbl_ldsc_1 = tk.Label(frame_inputs, text="N for Disease 1", font=("Roboto", round(user_height/70)))
    lbl_ldsc_2 = tk.Label(frame_inputs, text="N for Disease 2", font=("Roboto", round(user_height/70)))
    lbl_ldsc_3 = tk.Label(frame_inputs, text="Statistic (disease 1)", font=("Roboto", round(user_height/70)))
    lbl_ldsc_4 = tk.Label(frame_inputs, text="Statistic (disease 2)", font=("Roboto", round(user_height/70)))
    vcmd = root.register(integer_validate)
    ldsc_input_1 = tk.Entry(frame_inputs, textvariable=ldsc_n1, validate="all", validatecommand=(vcmd, "%P"))
    ldsc_input_2 = tk.Entry(frame_inputs, textvariable=ldsc_n2, validate="all", validatecommand=(vcmd, "%P"))
    ldsc_input_3 = tk.Entry(frame_inputs, textvariable=ldsc_stat1)
    ldsc_input_4 = tk.Entry(frame_inputs, textvariable=ldsc_stat2)

    result_label1 = tk.Label(frame_outputs, font=("Roboto", round(user_height/35), "bold"))
    result_label2 = tk.Label(frame_outputs)
    result_label3 = tk.Label(frame_outputs)

    # Calling the widgets
    # Title
    frame_title.pack()
    lbl_title.pack(anchor="center")

    # Content
    # Inputs
    frame_content.pack(anchor="center")
    frame_inputs.grid(row=0, column=0)

    frame_inputs.columnconfigure(0, weight=1, minsize=(user_width*0.25*0.50))
    frame_inputs.columnconfigure(1, weight=1, minsize=(user_width*0.25*0.50))


    lbl_inputs_0.grid(row=1,column=0, columnspan=2, sticky="w")
    upload_1.grid(row=2, column=0, columnspan=1, sticky="nsew")
    upload_2.grid(row=3, column=0, columnspan=1, sticky="nsew")
    lbl_inputs_1.grid(row=2, column=1, columnspan=2, sticky="w", padx=5)
    lbl_inputs_2.grid(row=3, column=1, columnspan=2, sticky="w", padx=5)
    frame_inputs.rowconfigure(4, minsize=30)
    lbl_inputs_3.grid(row=5,column=0, columnspan=2, sticky="w")
    radio1.grid(row=6,column=0, columnspan=2, sticky="w")
    radio2.grid(row=7,column=0, columnspan=2, sticky="w")
    radio4.grid(row=9,column=0, columnspan=2, sticky="w")
    frame_inputs.rowconfigure(10, minsize=30)
    lbl_inputs_4.grid(row=11,column=0, columnspan=2, sticky="w")
    frame_inputs.rowconfigure(15, minsize=20)
    go.grid(row=17,column=0, columnspan=2, sticky="w")

    lbl_ldsc_1.grid(row=12,column=0, columnspan=2, sticky="w", padx=5)
    lbl_ldsc_2.grid(row=13,column=0, columnspan=2, sticky="w", padx=5)
    lbl_ldsc_3.grid(row=14,column=0, columnspan=2, sticky="w", padx=5)
    lbl_ldsc_4.grid(row=15,column=0, columnspan=2, sticky="w", padx=5)
    ldsc_input_1.grid(row=12,column=1, columnspan=4, sticky="w")
    ldsc_input_2.grid(row=13,column=1, columnspan=4, sticky="w")
    ldsc_input_3.grid(row=14,column=1, columnspan=4, sticky="w")
    ldsc_input_4.grid(row=15,column=1, columnspan=4, sticky="w")
    show_test_variables()

    disease1_path = StringVar()
    disease2_path = StringVar()

    # Outputs
    frame_outputs.grid(row=0, column=1)
    lbl_outputs_1.pack(anchor="center")
    frame_outputs.columnconfigure(0, weight=1, minsize=(user_width*0.60))

    lbl_outputs_1.grid(row=1, column=0,  sticky="we")
    result_label1.grid(row=2,column=0, columnspan=10, sticky="we")
    result_label2.grid(row=3,column=0, columnspan=5, sticky="w", padx=5)
    result_label3.grid(row=3,column=0, columnspan=5, sticky="e", padx=5)
    # result_label.grid_remove()

    # Run the application
    root.mainloop()