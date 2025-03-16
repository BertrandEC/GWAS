# Importing module
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from pathlib import Path

# Defining Function
# Function to upload a file and display its content in the GUI
def upload_file(upload_label):
    file_path = filedialog.askopenfilename(filetypes=[("TSV files", "*.png")])
    if file_path:
        confirm = messagebox.askyesno("Confirm Upload", f"Are you sure you want to upload this file?\n{file_path}")
        if confirm:
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
            print(r)
            return ".../"+r
        else:
            return ".../"+q
    else:
        str(p)

def GO():
    pass

def check_files_uploaded():
    if lbl_inputs_1.cget("text") != "No file selected" and lbl_inputs_2.cget("text") != "No file selected":
        go.config(bg="green")


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
lbl_inputs_0 = tk.Label(frame_inputs, text="Upload GWAS Summary (.tsv.gz):", font=("Roboto", round(user_height/45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)
lbl_inputs_1 = tk.Label(frame_inputs, text="No file selected", font=("Roboto", round(user_height/70)))
lbl_inputs_2 = tk.Label(frame_inputs, text="No file selected", font=("Roboto", round(user_height/70)))
lbl_inputs_3 = tk.Label(frame_inputs, text="Select Statistical Method:", font=("Roboto", round(user_height/45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)
lbl_inputs_4 = tk.Label(frame_inputs, text="Select Correlation Test:", font=("Roboto", round(user_height/45), "bold"), fg="#333", bg="#B6D4F6", width=user_width)

# Defining button in frame_inputs
upload_1 = tk.Button(frame_inputs, text="Genetic Disease 1", font=("Roboto", round(user_height/70)), command=lambda: upload_file(lbl_inputs_1))
upload_2 = tk.Button(frame_inputs, text="Genetic Disease 2", font=("Roboto", round(user_height/70)), command=lambda: upload_file(lbl_inputs_2))
go = tk.Button(frame_inputs, text="GO!", font=("Roboto", round(user_height/60), "bold"), bg="#FF3B3B", fg="black", padx=20, pady=10, command=lambda: GO())

# Defining radio buttons in frame_inputs
stats_test = tk.StringVar(value="LR")
cor_method = tk.StringVar(value="P")
radio1 = tk.Radiobutton(frame_inputs, text="Linear Regression", font=("Roboto", round(user_height/70)), variable=stats_test, value="LR")
radio2 = tk.Radiobutton(frame_inputs, text="Linkage Disequilibrium Score Regression (LDSC)", font=("Roboto", round(user_height/70)), variable=stats_test, value="LDSC")
radio3 = tk.Radiobutton(frame_inputs, text="Pearson's", font=("Roboto", round(user_height/70)), variable=cor_method, value="P")
radio4 = tk.Radiobutton(frame_inputs, text="Spearman's", font=("Roboto", round(user_height/70)), variable=cor_method, value="S")


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


lbl_inputs_0.grid(row=1,column=0, columnspan=3, sticky="w")
upload_1.grid(row=2, column=0, columnspan=1, sticky="nsew")
upload_2.grid(row=3, column=0, columnspan=1, sticky="nsew")
lbl_inputs_1.grid(row=2, column=1, columnspan=2, sticky="w", padx=5)
lbl_inputs_2.grid(row=3, column=1, columnspan=2, sticky="w", padx=5)
frame_inputs.rowconfigure(4, minsize=30)
lbl_inputs_3.grid(row=5,column=0, columnspan=2, sticky="w")
radio1.grid(row=6,column=0, columnspan=4, sticky="w")
radio2.grid(row=7,column=0, columnspan=4, sticky="w")
frame_inputs.rowconfigure(8, minsize=30)
lbl_inputs_4.grid(row=9,column=0, columnspan=4, sticky="w")
radio3.grid(row=10,column=0, columnspan=4, sticky="w")
radio4.grid(row=11,column=0, columnspan=4, sticky="w")
frame_inputs.rowconfigure(12, minsize=20)
go.grid(row=13,column=0, columnspan=3, sticky="w")

# Outputs
frame_outputs.grid(row=0, column=1)

# Run the application
root.mainloop()