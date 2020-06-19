cd '\\carbon.research.sickkids.ca\rkafri\OPRETTA\Operetta Image Processing\User Interface\lib\NewVersion'
DirName='\\carbon.research.sickkids.ca\rkafri\OPRETTA\Operetta Processed OutPutFiles\Dataset_20170222_STRADa_SE_EdenRESULTS';
FileName='SegPar_SegmentationParameters';
load([DirName '\' FileName])
[ResultTable]=O_SegmentCells_v6(Data);
