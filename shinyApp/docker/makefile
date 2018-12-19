build:
	docker build -t hamid_tce -f Dockerfile .

bash:
	docker run -ti --rm -p 60040:3838 hamid_tce bash

dockerRun:
	docker run -p 60042:3838 hamid_tce

save:
	time docker save hamid_tce >  ~/tmp/docker_hamid_tce
	ls -lat ~/tmp/docker_hamid_tce

scp:
	time scp ~/tmp/docker_hamid_tce pshannon@trena:/local/tmp/
	ssh pshannon@trena ls -l /local/tmp/docker_hamid_tce
#	ssh pshannon@trena docker load < /local/tmp/dockerIGAP.$(ID)


