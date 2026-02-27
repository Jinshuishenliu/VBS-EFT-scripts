# Produce datacard
python3 makecard-boost-eft.py \
  --name eft_ZZ2l2nu \
  --input ./config/input_UL_2018-edirne-allDY-SR-aQGC.yaml \
  --era 2018 \
  --variable dilep_mt \
  --channel vbs-SR \
  --rebin=1 \
  --bins 100,300,400,500,600,700,800,900,1000

python3 makecard-boost-eft.py \
  --name eft_ZZ2l2nu \
  --input ./config/input_UL_2017-edirne-allDY-SR-aQGC.yaml \
  --era 2017 \
  --variable dilep_mt \
  --channel vbs-SR \
  --rebin=1 \
  --bins 100,300,400,500,600,700,800,900,1000

python3 makecard-boost-eft.py \
  --name eft_ZZ2l2nu \
  --input ./config/input_UL_2016-edirne-allDY-SR-aQGC.yaml \
  --era 2016 \
  --variable dilep_mt \
  --channel vbs-SR \
  --rebin=1 \
  --bins 100,300,400,500,600,700,800,900,1000

python3 makecard-boost-eft.py \
  --name eft_ZZ2l2nu \
  --input ./config/input_UL_2016APV-edirne-allDY-SR-aQGC.yaml \
  --era 2016APV \
  --variable dilep_mt \
  --channel vbs-SR \
  --rebin=1 \
  --bins 100,300,400,500,600,700,800,900,1000
