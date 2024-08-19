from waitress import serve
from app import server #"app" is the name of my Dash script I want to serve

serve(server)