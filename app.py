from flask import Flask, render_template, request

app = Flask(__name__)

@app.route('/', methods=['GET', 'POST'])
def index():
    selected_option = None
    if request.method == 'POST':
        selected_option = request.form.get('option')
    
    # Options for radio buttons
    options = [
        {'id': 'option1', 'value': 'Option 1', 'label': 'Option 1'},
        {'id': 'option2', 'value': 'Option 2', 'label': 'Option 2'},
        {'id': 'option3', 'value': 'Option 3', 'label': 'Option 3'}
    ]
    
    return render_template('index.html', options=options, selected_option=selected_option)

if __name__ == '__main__':
    app.run(debug=True)