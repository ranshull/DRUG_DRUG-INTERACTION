from flask import Flask, render_template, request, redirect, url_for
import joblib
from tensorflow.keras.models import load_model
import numpy as np
import cv2
import os
from flask import Flask, render_template, request, jsonify
import base64



app = Flask(__name__, static_folder='static')

#------------------------------------------------------------------
# Load the trained mode
model_path = 'models/2_stress_model.joblib'                                         
model = joblib.load(model_path)
#--------------------------------------------------------------------

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/audio-therapy')
def audio_therapy():
    return render_template('audioTherapy.html')

@app.route('/yoga-therapy')
def yoga_therapy():
    return render_template('yogatherapy.html')

@app.route('/PandE')
def prediction():
    return render_template('PandE.html')

@app.route('/contact')
def contact():
    return render_template('contact.html')

@app.route('/login')
def login():
    return render_template('login.html')

@app.route('/logout', methods=['GET'])
def logout():
    # Remove user session
    session.pop('username', None)
    return jsonify({'status': 'success', 'message': 'Logged out successfully!'}), 200

#------------------------------------------------------------

@app.route('/stress')
def Pstress():
    return render_template('StressPrediction.html')

@app.route('/submit_form', methods=['POST'])
def submit_form():
     # Retrieve form data
    rr = float(request.form['rr'])
    t = float(request.form['t'])
    lm = float(request.form['lm'])
    bo = float(request.form['bo'])
    rem = float(request.form['rem'])
    sr1 = float(request.form['sr1'])
    hr = float(request.form['hr'])
    
   # Create an input array for the model
    input_data = np.array([[rr, t, lm, bo, rem, sr1, hr]])

    # Make prediction
    prediction = model.predict(input_data)[0]  # This will be a numerical value (0-4)
    
    # Map the prediction to a stress level
    stress_levels = ["Low/Normal", "Medium Low", "Medium", "Medium High", "High"]
    stress_category = stress_levels[prediction]

       # Determine recommended therapy based on stress level
    if prediction >= 3:  # High or Medium High Stress
        therapy_link = url_for('yoga_therapy')  # Correctly link to yoga therapy page
        therapy_name = "Yoga"
    else:  # Low, Medium Low, or Medium Stress
        therapy_link = url_for('audio_therapy')  # Correctly link to audio therapy page
        therapy_name = "Audio"


   # Render the results page with the prediction and therapy link
    return render_template('2_result.html', prediction=prediction, stress_category=stress_category, therapy_link=therapy_link, therapy_name=therapy_name)

#-----------------------------------------------------------------------------------------------------------------------------------------------------

@app.route('/emotion')
def emotion():
    return render_template('emotiondetection.html')

# Load the trained emotion model
model_path = 'models\emotion_model.h5'
emotion_model = load_model(model_path)

# Emotion labels in the same order as your dataset
emotion_labels = ['Angry', 'Disgust', 'Fear', 'Happy', 'Neutral', 'Sad', 'Surprise']

# Add this feedback dictionary
emotion_feedback = {
    'Angry': 'Take a few deep breaths and try to relax.',
    'Disgust': 'Remember to focus on the positives around you.',
    'Fear': 'Stay calm; face your fears one step at a time.',
    'Happy': 'Keep smiling! Happiness is contagious.',
    'Neutral': 'All set; stay balanced and steady.',
    'Sad': 'Consider talking to a friend or doing something you enjoy.',
    'Surprise': 'Embrace the unexpected; life is full of surprises!'
}


@app.route('/predict_emotion', methods=['POST'])
def predict_emotion():
    if 'file' not in request.files:
        return jsonify({'error': 'No file uploaded'}), 400

    file = request.files['file']

    # Save and preprocess the uploaded file
    file_path = 'temp.jpg'
    file.save(file_path)

    img = cv2.imread(file_path, cv2.IMREAD_GRAYSCALE)
    img = cv2.resize(img, (48, 48))
    img = img.reshape(1, 48, 48, 1) / 255.0

    # Predict emotion
    predictions = emotion_model.predict(img)
    emotion_label = emotion_labels[np.argmax(predictions)]
    feedback_message = emotion_feedback[emotion_label]
    os.remove(file_path)  # Clean up temporary file
    return jsonify({
        'emotion': emotion_label,
        'feedback': feedback_message
    })


# @app.route('/predict_emotion', methods=['POST'])
# def predict_emotion():
#     # Parse the base64-encoded image from the POST request
#     data = request.get_json()
#     if 'image' not in data:
#         return jsonify({'error': 'No image data provided'}), 400

#     image_data = data['image']
#     image_data = image_data.split(",")[1]  # Remove the base64 header
#     img_bytes = base64.b64decode(image_data)

#     # Save and preprocess the image
#     file_path = 'temp.jpg'
#     with open(file_path, 'wb') as f:
#         f.write(img_bytes)

#     img = cv2.imread(file_path, cv2.IMREAD_GRAYSCALE)
#     img = cv2.resize(img, (48, 48))
#     img = img.reshape(1, 48, 48, 1) / 255.0

#     # Predict emotion
#     predictions = emotion_model.predict(img)
#     emotion_label = emotion_labels[np.argmax(predictions)]
#     feedback_message = emotion_feedback[emotion_label]
#     os.remove(file_path)  # Clean up temporary file
#     return jsonify({
#         'emotion': emotion_label,
#         'feedback': feedback_message
#     })

#---------------------------------------------------------------------------------------------------------------------------------

# from flask import Flask, render_template, request, jsonify, session, redirect, url_for
# from flask_sqlalchemy import SQLAlchemy
# import re  # For email validation


# # Initialize app and configure SQLite database
# app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'
# app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False

# db = SQLAlchemy(app)

# class User(db.Model):
#     id = db.Column(db.Integer, primary_key=True)
#     username = db.Column(db.String(80), unique=True, nullable=False)
#     email = db.Column(db.String(120), unique=True, nullable=False)
#     password = db.Column(db.String(120), nullable=False)

# app.secret_key = b'q\x82_\xe0d\xb2\xb7@y\x11&M\x13\x1d\xed\x8c\xed\xab\x9e\x81\x1d\xe8m\xae'


# @app.route('/signup', methods=['POST'])
# def signup():
#     data = request.get_json()
#     print("Received signup data:", data)

#     username = data.get('signupUsername')
#     email = data.get('signupEmail')
#     password = data.get('signupPassword')

#     if not username or not email or not password:
#         return jsonify({'status': 'error', 'message': 'Please fill in all fields!'}), 400

#     # Check if user already exists
#     existing_user = User.query.filter_by(email=email).first()
#     if existing_user:
#         return jsonify({'status': 'error', 'message': 'Email already registered!'}), 409

#     new_user = User(username=username, email=email, password=password)
#     db.session.add(new_user)
#     db.session.commit()

#     return jsonify({'status': 'success', 'message': 'Signup successful!'}), 201


# @app.route('/login_user', methods=['POST'])
# def login_user():
#     data = request.get_json()
#     email = data.get('loginEmail')
#     password = data.get('loginPassword')

#     user = User.query.filter_by(email=email, password=password).first()
#     if user:
#         return jsonify({'status': 'success', 'message': 'Login successful!'}), 200
#     else:
#         return jsonify({'status': 'error', 'message': 'Invalid credentials!'}), 401

# from flask import Flask, render_template, request, jsonify
# from flask_pymongo import PyMongo

# app = Flask(__name__)

# # Configure MongoDB connection
# app.config["MONGO_URI"] = "mongodb://localhost:27017/mydatabase"  # Replace with your MongoDB URI
# mongo = PyMongo(app)

# # Access the MongoDB collection
# users_collection = mongo.db.users



# # Example to add a user
# @app.route('/signup', methods=['POST'])
# def signup():
#     data = request.get_json()
#     username = data.get('username')
#     email = data.get('email')
#     password = data.get('password')
    
#     # Check if user already exists
#     if users_collection.find_one({"email": email}):
#         return jsonify({"status": "error", "message": "User already exists!"}), 409

#     # Insert the new user
#     users_collection.insert_one({"username": username, "email": email, "password": password})
#     return jsonify({"status": "success", "message": "Signup successful!"}), 201

# # Example to retrieve user for login
# @app.route('/login', methods=['POST'])
# def login():
#     data = request.get_json()
#     email = data.get('email')
#     password = data.get('password')
    
#     # Check credentials
#     user = users_collection.find_one({"email": email, "password": password})
#     if user:
#         return jsonify({"status": "success", "message": "Login successful!"}), 200
#     else:
#         return jsonify({"status": "error", "message": "Invalid credentials!"}), 401

# if __name__ == "__main__":
#     app.run(debug=True)


# ----------------------------------------------------------

# # updated code for db
# from flask import Flask, render_template, request, jsonify, session, redirect, url_for
# from flask_sqlalchemy import SQLAlchemy
# import re  # For email validation
# from werkzeug.security import generate_password_hash, check_password_hash

# # Initialize app
# app = Flask(__name__, static_folder='static')

# # Configurations
# app.config['SQLALCHEMY_DATABASE_URI'] = 'sqlite:///users.db'
# app.config['SQLALCHEMY_TRACK_MODIFICATIONS'] = False
# app.secret_key = b'q\x82_\xe0d\xb2\xb7@y\x11&M\x13\x1d\xed\x8c\xed\xab\x9e\x81\x1d\xe8m\xae'

# # Initialize the database
# db = SQLAlchemy(app)

# # Define User Model
# class User(db.Model):
#     id = db.Column(db.Integer, primary_key=True)
#     username = db.Column(db.String(80), unique=True, nullable=False)
#     email = db.Column(db.String(120), unique=True, nullable=False)
#     password = db.Column(db.String(120), nullable=False)

# # Home Page
# @app.route('/')
# def home():
#     username = session.get('username')  # Check if the user is logged in
#     return render_template('index.html', username=username)

# # Sign-Up Route
# @app.route('/signup', methods=['POST'])
# def signup():
#     data = request.get_json()
#     username = data.get('signupUsername')
#     email = data.get('signupEmail')
#     password = data.get('signupPassword')

#     # Validate email format
#     if not re.match(r"[^@]+@[^@]+\.[^@]+", email):
#         return jsonify({'status': 'error', 'message': 'Invalid email format!'}), 400

#     # Check if email already exists
#     existing_user = User.query.filter_by(email=email).first()
#     if existing_user:
#         return jsonify({'status': 'error', 'message': 'Email already registered!'}), 409

#     # Hash the password
#     hashed_password = generate_password_hash(password)

#     # Create and save the new user
#     new_user = User(username=username, email=email, password=hashed_password)
#     db.session.add(new_user)
#     db.session.commit()
#     return jsonify({'status': 'success', 'message': 'Signup successful!'}), 201

# # Login Route
# @app.route('/login_user', methods=['POST'])
# def login_user():
#     data = request.get_json()
#     email = data.get('loginEmail')
#     password = data.get('loginPassword')

#     # Check if user exists
#     user = User.query.filter_by(email=email).first()
#     if user and check_password_hash(user.password, password):
#         session['username'] = user.username  # Store username in session
#         return jsonify({'status': 'success', 'message': 'Login successful!'}), 200
#     else:
#         return jsonify({'status': 'error', 'message': 'Invalid credentials!'}), 401

# # Logout Route
# @app.route('/logout', methods=['GET'])
# def logout():
#     session.pop('username', None)  # Remove username from session
#     return redirect(url_for('home'))



if __name__ == '__main__':
    app.run(debug=True)

# python app.py
