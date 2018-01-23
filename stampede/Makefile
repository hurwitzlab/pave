APP = imicrobe-mash-0.0.1

test: clean
	sbatch test.sh

clean:
	find . \( -name \*.conf -o -name \*.out -o -name \*.log -o -name \*.params -o -name launcher-\* \) -exec rm {} \;

upload:
	files-upload -F $(WORK)/$(APP) kyclark/applications

update:
	apps-addupdate -F app.json
