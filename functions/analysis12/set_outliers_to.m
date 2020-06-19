function data = set_outliers_to(data, lower_prctile, higher_prctile, new_value)
  lo_prctile = prctile(data,lower_prctile);
  hi_prctile = prctile(data,higher_prctile);
  data(data < lo_prctile | data > hi_prctile) = new_value;
end