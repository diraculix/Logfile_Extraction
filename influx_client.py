import os, platform, psutil, time
from datetime import datetime
from influxdb_client import InfluxDBClient, Point, WritePrecision
from influxdb_client.client.write_api import SYNCHRONOUS

# access config
token = "vHLCAHK5mURBcgOJane0P6GsOyBkaDg1H9cHTGwgEwWlRwS7UuFitjIewPhf_o9AZVwinuzHyDD6c9bL6HN2Ug=="
org = "OncoRay"
bucket = "GTR2"
url = "http://g40invivodossrv:8086"

# get host name
if platform.system() == "Windows":
    host = platform.uname().node
else:
    host = os.uname()[1]  # does not work on windows

# monitor loop
print(f'[{time.ctime()}] Establishing connection to {url}/{org}/{bucket} as {host}..')
with InfluxDBClient(url=url, token=token, org=org) as client:
    write_api = client.write_api(write_options=SYNCHRONOUS)
    while True:
        print(f'[{time.ctime()}] > Influx: Connected to database, monitor running..', end='\r')
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
            # time.sleep(5)
        
        except KeyboardInterrupt:
            print(f'\n[{time.ctime()}] > Influx: Monitor on host {host} stopped, disconnecting client..')
            break
    
    client.close()
