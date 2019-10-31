StateMWGwasPeakHyper.jar: manifest.txt
	javac -Xlint ./StateMWGwasPeakHyper.java
		jar cvfm StateMWGwasPeakHyper.jar manifest.txt StateMWGwasPeakHyper*.class

manifest.txt: 
		echo 'Main-Class: StateMWGwasPeakHyper' > manifest.txt

clean:
	rm manifest.txt StateMWGwasPeakHyper.jar
