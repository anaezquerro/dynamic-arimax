build:
	docker build . -t dynamic-arimax

run: 
	docker run --rm -d -it --name dyntest --mount type=bind,source=C:\Users\ana\Documents\dynamic-arimax\,target=/app dynamic-arimax
	docker attach dyntest

stop:
	docker stop dyntest
	docker rm dyntest

rm:
	docker image rm dynamic-arimax