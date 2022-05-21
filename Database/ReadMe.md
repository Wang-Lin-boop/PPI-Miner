Database
====

1. Decompress the zip file.

```
unzip 3DMotif-Dock-pds-Library.zip
```

2. Generate index file.

```
for pds in `ls 3DMotif-Dock-pds-Library`;do
  echo "${PWD}/3DMotif-Dock-pds-Library/${pds}" >> 3DMotif-Dock-pds-Library.list
done
```
