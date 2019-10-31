ChromImpute.jar: manifest.txt
	javac -Xlint ./ernst/ChromImpute/*.java
	jar cvfm ChromImpute.jar manifest.txt ernst/ChromImpute/*.class

manifest.txt: 
	echo 'Main-Class: ernst.ChromImpute.ChromImpute' > manifest.txt

clean:
	rm manifest.txt ChromImpute.jar
