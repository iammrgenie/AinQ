#For the Team Leader Drone
TL-AinQ: TL-ainq.c
	gcc -O2 -std=c99 TL-ainq.c core.a -o TL-AinQ

#For the Edge Devices
#ED-AinQ: ED-ainq.c
#	gcc -O2 -std=c99 ED-ainq.c core.a -o ED-AinQ

.PHONY: clean

clean:
	rm -rf TL-AinQ