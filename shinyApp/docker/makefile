ID=19dec2018

build:
	docker build -t hamid_tce$(ID) -f Dockerfile .

bash:
	docker run -ti --rm -p 60040:3838 hamid_tce$(ID) bash

dockerRun:
	docker run -p 60042:3838 hamid_tce$(D)

save:
	time docker save hamid_tce$(ID) >  ~/tmp/docker_hamid_tce$(ID)
	ls -lat ~/tmp/docker_hamid_tce$(ID)

scp:
	time scp ~/tmp/docker_hamid_tce$(ID) pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/docker_hamid_tce
#	ssh pshannon@trena docker load < /local/tmp/dockerIGAP.$(ID)


