# Importing module
import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox
from pathlib import Path




# Defining Function
# Function to upload a file and display its content in the GUI
def upload_file(upload_label):
    file_path = filedialog.askopenfilename(filetypes=[("GWAS summary statistics", "*.jpg")])
    if file_path:
        confirm = messagebox.askyesno("Confirm Upload", f"Are you sure you want to upload this file?\n{file_path}")
        if confirm:
            upload_label.config(text=f"{shorten_path(file_path)} uploaded")
        else:
            upload_label.config(text="Upload Cancelled")
            upload_label.config(text="No file selected")

def shorten_path(full_path):
    p = Path(full_path)
    return ".../"+"/".join(p.parts[-2:]) if len(p.parts) > 2 else str(p)

def GO():
    pass

# Setting up the window
# Defining an instance window called "root" into Tk class
root = tk.Tk()
root.title("CLAS: Cross-Linking Association Studies")
root.geometry("700x700")
# Getting user's window information
user_width = root.winfo_screenwidth()
user_height = root.winfo_screenheight()



# Defining Widgets
# Create frame for title: frame_title
frame_title = tk.Frame(root, bg="#000080", width=user_width, height=user_height/12)
frame_title.pack_propagate(False)  # Prevent the frame from shrinking to fit its content
# Create frame for content: frame_content
frame_content = tk.Frame(root, width=user_width, height=(user_height-(user_height/12)))
frame_content.grid_propagate(False)
#Create frame for inputs: frame_inputs
frame_inputs = tk.Frame(frame_content, width=(user_width*0.25), height=(user_height-(user_height/12)))
frame_inputs.grid_propagate(False)
#Create frame for outputs: frame_outputs
frame_outputs = tk.Frame(frame_content, width=(user_width*0.75), height=(user_height-(user_height/12)))
frame_outputs.grid_propagate(False)


# Define title in the frame_title
lbl_title = tk.Label(frame_title, text="CLAS: Cross-Linking Association Studies",
                     font=("Arial", round(user_height/25)), fg="white", bg="#000080", anchor="center")

# Defining label in frame_inputs
lbl_inputs_0 = tk.Label(frame_inputs, text="Upload GWAS results in .tsv.gz:", font=("Arial", round(user_height/50)))
lbl_inputs_1 = tk.Label(frame_inputs, text="No file selected", font=("Arial", round(user_height/100)))
lbl_inputs_2 = tk.Label(frame_inputs, text="No file selected", font=("Arial", round(user_height/100)))
lbl_inputs_3 = tk.Label(frame_inputs, text="Select statistical method:", font=("Arial", round(user_height/50)))
lbl_inputs_4 = tk.Label(frame_inputs, text="Select correlation test:", font=("Arial", round(user_height/50)))
# Defining button in frame_inputs
upload_1 = tk.Button(frame_inputs, text="Genetic Disease 1", font=("Arial", round(user_height/100)), command=lambda: upload_file(lbl_inputs_1))
upload_2 = tk.Button(frame_inputs, text="Genetic Disease 2", font=("Arial", round(user_height/100)), command=lambda: upload_file(lbl_inputs_2))
go = tk.Button(frame_inputs, text="GO!", bg="green", padx=30, pady=10, command=lambda: GO())
# Defining radio buttons in frame_inputs
stats_test = tk.StringVar(value="")
cor_method = tk.StringVar(value="")
radio1 = tk.Radiobutton(frame_inputs, text="Linear Regression", font=("Arial", round(user_height/100)) , variable=stats_test, value="LR")
radio2 = tk.Radiobutton(frame_inputs, text="LDSC", font=("Arial", round(user_height/100)), variable=stats_test, value="LDSC")
radio3 = tk.Radiobutton(frame_inputs, text="Pearson's", font=("Arial", round(user_height/100)), variable=cor_method, value="P")
radio4 = tk.Radiobutton(frame_inputs, text="Spearman's", font=("Arial", round(user_height/100)), variable=cor_method, value="S")


# Calling the widgets
# Title
frame_title.pack()
lbl_title.pack(anchor="center")

# Content
# Inputs
frame_content.pack(anchor="center")
frame_inputs.grid(row=0, column=0)

frame_inputs.columnconfigure(0, weight=1, minsize=(user_width*0.25*0.40))
frame_inputs.columnconfigure(1, weight=1, minsize=(user_width*0.25*0.60))

lbl_inputs_0.grid(row=0,column=0, columnspan=2, sticky="w")
upload_1.grid(row=1, column=0, sticky="nsew")
upload_2.grid(row=2, column=0, sticky="nsew")
lbl_inputs_1.grid(row=1, column=1, sticky="w")
lbl_inputs_2.grid(row=2, column=1, sticky="w")
frame_inputs.rowconfigure(3, minsize=30)
lbl_inputs_3.grid(row=4,column=0, columnspan=2, sticky="w")
radio1.grid(row=5,column=0, sticky="w")
radio2.grid(row=6,column=0, sticky="w")
frame_inputs.rowconfigure(7, minsize=30)
lbl_inputs_4.grid(row=8,column=0, columnspan=2, sticky="w")
radio3.grid(row=9,column=0, sticky="w")
radio4.grid(row=10,column=0, sticky="w")
frame_inputs.rowconfigure(11, minsize=20)
go.grid(row=12,column=0, columnspan=3, sticky="w")

# Outputs
frame_outputs.grid(row=0, column=1)

# Run the application
root.mainloop()