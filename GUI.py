import tkinter as tk
from tkinter import ttk
from tkinter import filedialog
from tkinter import messagebox # Import the messagebox module



# Defining window called "root"
root = tk.Tk()
root.title("CLAS: Cross-Linking Association Studies")
root.geometry("500x500")

# Getting user's window information
user_width = root.winfo_screenwidth()
user_height = root.winfo_screenheight()



# Create a frame for title: frame_title
frame_title = tk.Frame(root, bg="#000080", width=user_width, height=user_height/12)
frame_title.pack_propagate(False)  # Prevent the frame from shrinking to fit its content
frame_title.pack()  # Pack with padding around the frame

# Create a title in the frame_title
lbl_title = tk.Label(frame_title, text="CLAS: Cross-Linking Association Studies",
                     font=("Arial", round(user_height/25)), fg="white", bg="#000080")
lbl_title.pack(anchor="center")



# Generating frame for inputs: frame_inputs
frame_inputs = tk.Frame(root, width=user_width, height=(user_height-(user_height/12)))
frame_inputs.pack_propagate(False)
frame_inputs.pack()

# Create a label to display the file path
label_1 = tk.Label(frame_inputs, text="No file selected")
label_1.grid(row=0, column=0, pady=10, padx=50, sticky="nsew")

label_2 = tk.Label(frame_inputs, text="No file selected")
label_2.grid(row=0, column=1, pady=10, padx=50, sticky="nsew")



# Function to upload a file and display its content in the GUI
def upload_file(upload_label):
    file_path = filedialog.askopenfilename(filetypes=[("GWAS summary statistics", "*.tsv.gz")])
    if file_path:
        confirm = messagebox.askyesno("Confirm Upload", f"Are you sure you want to upload this file?\n{file_path}")
        if confirm:
            upload_label.config(text=f"File selected: {file_path} uploaded")
        else:
            upload_label.config(text="Upload Cancelled")
            upload_label.config(text="No file selected")



# Create a button to trigger the file upload
upload_1 = tk.Button(frame_inputs, text="Upload Dataframe for Genetic Disease 1", command=lambda: upload_file(label_1))
upload_1.grid(row=1, column=0, pady=10, padx=50, sticky="nsew")

upload_2 = tk.Button(frame_inputs, text="Upload Dataframe for Genetic Disease 2", command=lambda: upload_file(label_2))
upload_2.grid(row=1, column=1, pady=10, padx=50, sticky="nsew")

# Run the application
root.mainloop()