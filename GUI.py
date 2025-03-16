# Importing module
import tkinter as tk
from tkinter import filedialog, StringVar
from tkinter import messagebox
from pathlib import Path
import pandas as pd
import hashlib
from threading import Thread

from gwas_analysis.plots import plot_joint_manhattan, plot_scatter_with_regression
from gwas_analysis.stats import tidy_summary_stats

# Defining Function
# Function to upload a file and display its content in the GUI
def upload_file(upload_label):
    file_path = filedialog.askopenfilename(filetypes=[("TSV files", "*.tsv")])
    if file_path:
        # confirm = messagebox.askyesno("Confirm Upload", f"Are you sure you want to upload this file?\n{file_path}")
        confirm = True
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
    with open(disease1_path.get(), 'rb', buffering=0) as f:
        hash1 = hashlib.file_digest(f, 'sha1').hexdigest()
    with open(disease2_path.get(), 'rb', buffering=0) as f:
        hash2 = hashlib.file_digest(f, 'sha1').hexdigest()

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
    match stats_test.get():
        case "LR":
            filtered_df = merged_df[(merged_df['p_value_diabetes'] < 1e-5) | (merged_df['p_value_arthritis'] < 1e-5)]
            if cor_method.get() == "B":
                plot_scatter_with_regression(filtered_df, col1='beta_disease1', col2='beta_disease2')
            elif cor_method.get() == "P":
                plot_scatter_with_regression(filtered_df, col1='neg_log_p_value_disease1', col2='neg_log_p_value_disease2')
        case "LDSC": pass
        case "M":
            plot_joint_manhattan(merged_df, 'neg_log_p_value_disease1', 'neg_log_p_value_disease2')

    root.after(0, display_results)

def display_results():
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
    method1.grid_remove()
    method2.grid_remove()
    lbl_ldsc_1.grid_forget()
    lbl_ldsc_2.grid_forget()
    ldsc_input_1.grid_forget()
    ldsc_input_2.grid_forget()
    match stats_test.get():
        case "LR":
            method1.grid()
            method2.grid()
        case "LDSC":
            lbl_ldsc_1.grid(row=11,column=0, columnspan=2, sticky="w", padx=5)
            lbl_ldsc_2.grid(row=12,column=0, columnspan=2, sticky="w", padx=5)
            ldsc_input_1.grid(row=11,column=1, columnspan=4, sticky="w")
            ldsc_input_2.grid(row=12,column=1, columnspan=4, sticky="w")
        case "M":
            pass

if __name__ == "__main__":

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
    lbl_inputs_3 = tk.Label(frame_inputs, text="Select Statistical Method:", font=("Roboto", round(user_height/45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)
    lbl_inputs_4 = tk.Label(frame_inputs, text="Select Options:", font=("Roboto", round(user_height/45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)
    # lbl_inputs_5 = tk.Label(frame_inputs, text="Results:", font=("Roboto", round(user_height/45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)

    # Defining button in frame_inputs
    upload_1 = tk.Button(frame_inputs, text="Genetic Disease 1", font=("Roboto", round(user_height/70)), command=lambda: upload_file(lbl_inputs_1))
    upload_2 = tk.Button(frame_inputs, text="Genetic Disease 2", font=("Roboto", round(user_height/70)), command=lambda: upload_file(lbl_inputs_2))
    go = tk.Button(frame_inputs, text="GO!", font=("Roboto", round(user_height/60), "bold"), bg="#FF3B3B", fg="black", padx=20, pady=10, command=lambda: GO())

    # Defining radio buttons in frame_inputs
    stats_test = tk.StringVar(value="LR")
    radio1 = tk.Radiobutton(frame_inputs, text="Linear Regression", font=("Roboto", round(user_height/70)), variable=stats_test, value="LR", command=show_test_variables)
    radio2 = tk.Radiobutton(frame_inputs, text="Linkage Disequilibrium Score Regression (LDSC)", font=("Roboto", round(user_height/70)), variable=stats_test, value="LDSC", command=show_test_variables)
    radio3 = tk.Radiobutton(frame_inputs, text="Generate Manhattan Plot", font=("Roboto", round(user_height/70)), variable=stats_test, value="M", command=show_test_variables)

    cor_method = tk.StringVar(value="P")
    method1 = tk.Radiobutton(frame_inputs, text="Beta values", font=("Roboto", round(user_height / 70)), variable=cor_method, value="B")
    method2 = tk.Radiobutton(frame_inputs, text="P values", font=("Roboto", round(user_height / 70)), variable=cor_method, value="P")

    ldsc_n1 = tk.IntVar(value=0)
    ldsc_n2 = tk.IntVar(value=0)
    lbl_ldsc_1 = tk.Label(frame_inputs, text="N for Disease 1", font=("Roboto", round(user_height/70)))
    lbl_ldsc_2 = tk.Label(frame_inputs, text="N for Disease 2", font=("Roboto", round(user_height/70)))
    vcmd = root.register(integer_validate)
    ldsc_input_1 = tk.Entry(frame_inputs, textvariable=ldsc_n1, validate="all", validatecommand=(vcmd, "%P"))
    ldsc_input_2 = tk.Entry(frame_inputs, textvariable=ldsc_n2, validate="all", validatecommand=(vcmd, "%P"))

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
    # frame_inputs.columnconfigure(2, weight=1, minsize=(user_width*0.50))


    lbl_inputs_0.grid(row=1,column=0, columnspan=2, sticky="w")
    upload_1.grid(row=2, column=0, columnspan=1, sticky="nsew")
    upload_2.grid(row=3, column=0, columnspan=1, sticky="nsew")
    lbl_inputs_1.grid(row=2, column=1, columnspan=2, sticky="w", padx=5)
    lbl_inputs_2.grid(row=3, column=1, columnspan=2, sticky="w", padx=5)
    frame_inputs.rowconfigure(4, minsize=30)
    lbl_inputs_3.grid(row=5,column=0, columnspan=2, sticky="w")
    radio1.grid(row=6,column=0, columnspan=2, sticky="w")
    radio2.grid(row=7,column=0, columnspan=2, sticky="w")
    radio3.grid(row=8,column=0, columnspan=2, sticky="w")
    frame_inputs.rowconfigure(9, minsize=30)
    lbl_inputs_4.grid(row=10,column=0, columnspan=2, sticky="w")
    method1.grid(row=11, column=0, columnspan=2, sticky="w")
    method2.grid(row=12, column=0, columnspan=2, sticky="w")
    # lbl_inputs_5.grid(row=1, column=2, columnspan=2, sticky="w")
    frame_inputs.rowconfigure(13, minsize=20)
    go.grid(row=14,column=0, columnspan=2, sticky="w")

    disease1_path = StringVar()
    disease2_path = StringVar()

    # Outputs
    frame_outputs.grid(row=0, column=1)

    # Run the application
    root.mainloop()