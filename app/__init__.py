#!/usr/bin/python
# __init__.py
# This is the file that initializes the "Flask" object in python. This script
# basically initializes everything that the app needs to run properly.

from flask import Flask

app = Flask(__name__)
app.config.from_object('config')

from app import views
