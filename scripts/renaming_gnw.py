for file in *_dream4_timeseries.tsv;
  do mv "${file}" "${file/_dream4/}";
done
