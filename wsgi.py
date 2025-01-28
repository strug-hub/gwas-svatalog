from waitress import serve
from app import server 

serve(server, 
      listen = '127.0.0.1:1234', 
      threads = 10)