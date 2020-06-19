function fun(ResultTable, Data)
  stats = table();
  stats.num_timepoints= length(unique(ResultTable.Time));
  stats.num_rows_used = length(unique(ResultTable.Row));
  stats.num_columns_used = length(unique(ResultTable.Column));
  stats.num_fields_used = length(unique(ResultTable.Field));
  stats.channels = Data.O.Channels;
  stats.plate_type = Data.O.PlateType;
  stats.num_images = size(Data.O.ImageIDs,1);
  stats.num_cells = height(ResultTable);
  stats.magnification = Data.O.SegmentationParameters.magnification;
  stats.dataset_name = Data.O.SegmentationParameters.DataSetName;
  stats.raw_images_location = Data.O.SegmentationParameters.DirName;
  sample_image = imread([stats.raw_images_location '\' char(Data.O.ImageIDs.FileName(1))]);
  stats.image_resolution = size(sample_image);
  stats
  writetable(stats,'basic_stats.xls');
end