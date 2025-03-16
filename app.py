from flask import Flask , render_template
from Form.form import form
from SIMPLE.simpleChecker import simpleChecker 

app =  Flask(__name__)

app.register_blueprint(form, url_prefix='/form')
app.register_blueprint(simpleChecker , url_prefix='/simpleChecker')

@app.route("/")
def test():
    return render_template("main.html")

if __name__ == '__main__':
    app.run(debug=True)
