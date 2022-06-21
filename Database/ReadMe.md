Database
====

1. Download and decompress the zip file.

```
wget https://bailab.siais.shanghaitech.edu.cn/service/3DMotif-Dock-pds-Library.zip
unzip 3DMotif-Dock-pds-Library.zip
```

2. Generate index file.

```
for pds in `ls 3DMotif-Dock-pds-Library`;do
  echo "${PWD}/3DMotif-Dock-pds-Library/${pds}" >> 3DMotif-Dock-pds-Library.list
done
```
