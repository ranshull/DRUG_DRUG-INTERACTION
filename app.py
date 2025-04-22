from flask import Flask, render_template, request, jsonify  # ✅ Added request
from flask import Flask, render_template, request, redirect, flash
from pymongo import MongoClient
import re
from werkzeug.security import generate_password_hash, check_password_hash

from Form.form import form
from SIMPLE.simpleChecker import simpleChecker

app = Flask(__name__)
app.secret_key = "a_random_and_unique_secret_key_1234"




# from pymongo.mongo_client import MongoClient
# from pymongo.server_api import ServerApi

# uri = "mongodb+srv://ANSHUL-MEDSAFE:<Geetaanuj04>@cluster0.muostsi.mongodb.net/?retryWrites=true&w=majority&appName=Cluster0"

# # Create a new client and connect to the server
# client = MongoClient(uri, server_api=ServerApi('1'))

# # Send a ping to confirm a successful connection
# try:
#     client.admin.command('ping')
#     print("Pinged your deployment. You successfully connected to MongoDB!")
# except Exception as e:
#     print(e)

# MongoDB connection (local)
client = MongoClient("mongodb+srv://ANSHUL-MEDSAFE:Yl9B8Rl5QhXdK1D7@cluster0.muostsi.mongodb.net/?retryWrites=true&w=majority&appName=Cluster0")
db = client.user_auth
users = db.users

app.register_blueprint(form, url_prefix='/form')
app.register_blueprint(simpleChecker, url_prefix='/simpleChecker')

# Email regex
email_regex = re.compile(r'^\S+@\S+\.\S+$')

@app.route("/")
def test():
    return render_template("main.html")

# @app.route("/login")
# def login():
#     return render_template("login.html")

@app.route('/signup', methods=['GET', 'POST'])
def signup():
    if request.method == 'POST':
        name = request.form['name'].strip()
        email = request.form['email'].strip().lower()
        password = request.form['password']
        confirm_password = request.form['confirm_password']

        # Backend validations
        if not name or not email or not password or not confirm_password:
            flash("Please fill all fields")
        elif not email_regex.match(email):
            flash("Invalid email format")
        elif len(password) < 6:
            flash("Password must be at least 6 characters long")
        elif password != confirm_password:
            flash("Passwords do not match")
        elif users.find_one({"email": email}):
            flash("Email already registered")
        else:
            hashed_pw = generate_password_hash(password)
            users.insert_one({
                "name": name,
                "email": email,
                "password": hashed_pw
            })
            flash("Sign up successful! Please login.")
            return redirect('/login')

    return render_template('signup.html')

@app.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        email = request.form['email'].strip().lower()
        password = request.form['password']

        user = users.find_one({"email": email})
        if user and check_password_hash(user['password'], password):
            flash("Login successful!")
            return redirect('/')  # redirect to dashboard if needed
        else:
            flash("Invalid credentials")

    return render_template('login.html')

if __name__ == "__main__":
    app.run(debug=True)

