function fun(ResultTable, Data)
  channels = Data.O.Channels;
  raw_images_location = Data.O.SegmentationParameters.DirName;

  figure('Position',[251,    91,   983,   708])
  for i=1:length(channels)
    img = imread([raw_images_location '\' char(Data.O.ImageIDs.FileName(i))]);
    subplot(2,2,i)
    imshow(img, [prctile(img(:), 3), prctile(img(:), 99.9)])
    title(channels{i})
  end

  % Save the plot to disk
  filename = sprintf('plots\\sample_images.png');
  export_fig(filename, '-m4')
end