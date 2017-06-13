# PETPVC mask-related tools

## Relabeling an image

The program ```pvc_relabel``` can be used to transform a 3-D parcellation from region definition to another. The mapping is controlled by a CSV file which is passed as one of the inputs.

### Usage
```
    pvc_relabel -i <INPUT> -o <OUTPUT> --parc <CSVFILE> --type <COLID>
```
where ```<INPUT>``` is the 3-D image of parcellated regions, ```<OUTPUT>``` is the newly labeled image, ```<CSVFILE>``` is a comma separated file (see below) that maps from one labelling scheme to another and ```<COLID>``` refers to the desired mapping contained in a column of the ```<CSVFILE>```.

### CSV file definition

An example mapping from FreeSurfer's ```aparc+aseg``` is provided in ```parc/FS.csv```. ```pvc_relabel``` always expects the first row to comprise the column headers:

```
REGION,SOURCE,XXX
```

The ```REGION``` column contains a text description of the region (optional) and ```SOURCE``` holds the numeric ID of the original label (mandatory). The header ```XXX``` can be any desired name for the new mapping scheme. In ```FS.csv``` this is called ```DST```. The column ***must*** contain the new numeric ID to be applied. It is this name that must be passed as the ```<COLID>``` when executing ```pvc_relabel```. Multiple columns can be added to the CSV file for different applications:

```
REGION,SOURCE,MYSCHEME1,MYSCHEME2
```
Calling ```pvc_relabel``` with ```--type MYSCHEME1``` will produce an image with the label scheme defined by the ```MYSCHEME1``` and ```--type MYSCHEME2``` will use the second column.

Note that  ```FS.csv``` is an example and that region definitions should be created and validated for each application.

## Producing a 4-D mask file from 3-D labels

Most of the methods provided in the toolbox expect a 4-D mask image as input. The application ```pvc_make4d``` will create a 4-D image from a 3-D label image. Each volume in the 4-D file represents a single region. The set of labels do not need to be continuous. ```pvc_make4d``` will write the volumes in ascending order. 

### Usage
```
	pvc_make4d -i <3DMASK> -o <4DMASK>
```
where ```<3DMASK>``` is a single 3-D volume of labels and ```<4DMASK>``` is the output 4-D file.