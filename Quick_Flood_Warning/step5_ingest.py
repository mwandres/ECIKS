import pysftp
import os
rootDir='D:\\ECIKS\\Quick_Flood_Warning\\results\\'
Hostname = "192.168.53.68"
Username = "anujftp"
Password = "Simple10"
cnopts = pysftp.CnOpts()
cnopts.hostkeys = None

def ingest_products():
	for file in os.listdir(rootDir):
		with pysftp.Connection(host=Hostname, username=Username, password=Password,cnopts=cnopts) as sftp:
			print("Connection successfully established ... ")
			sftp.cwd('/home/anujftp')
			sftp.execute('rm -rf /home/anujftp/'+file)
			sftp.put(rootDir+file, '/home/anujftp/'+file)
	return True

