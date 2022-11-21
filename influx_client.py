import psutil
import time
from datetime import datetime

from influxdb_client import InfluxDBClient, Point, WritePrecision
from influxdb_client.client.write_api import SYNCHRONOUS

# You can generate an API token from the "API Tokens Tab" in the UI
token = "vHLCAHK5mURBcgOJane0P6GsOyBkaDg1H9cHTGwgEwWlRwS7UuFitjIewPhf_o9AZVwinuzHyDD6c9bL6HN2Ug=="
org = "OncoRay"
bucket = "GTR2"
url = "http://g40invivodossrv:8086"

print(f'[{datetime.utcnow()}] Establishing connection to {url}/{org}/{bucket} ..')
with InfluxDBClient(url=url, token=token, org=org) as client:
    host = "PC43doktoranden"
    write_api = client.write_api(write_options=SYNCHRONOUS)
    print(f'[{datetime.utcnow()}] > Influx: Connected to database, monitor running..', end='\r')

    while True:
        try:
            point = Point("cpu") \
                .tag("host", host) \
                .field("used_percent", psutil.cpu_percent(5)) \
                .time(datetime.utcnow(), WritePrecision.NS)

            write_api.write(bucket, org, point)

            point = Point("mem") \
                .tag("host", host) \
                .field("used_percent", psutil.virtual_memory()[2]) \
                .time(datetime.utcnow(), WritePrecision.NS)

            write_api.write(bucket, org, point)
            time.sleep(5)
        
        except KeyboardInterrupt:
            print(f'[{datetime.utcnow()}] > Influx: Monitor on host {host} stopped, disconnecting client..')
            break
    
    client.close()
