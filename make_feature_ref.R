##Make feature_ref.csv for input of cellanger multi for CITE-SEQ data (RNA+ADT)

library(readxl)
library(dplyr)
library(readr)

# Input Excel with columns: DNA_ID, Description, Barcode
#Use totalseqC.xlsx as input
xlsx <- "/net/beegfs/users/P089417/MDS_CRfiles/totalseqC.xlsx"

base <- read_excel(xlsx) %>%
  dplyr::select(DNA_ID, Description, Barcode) %>%
  dplyr::filter(!is.na(Barcode))

# Add your three spike-ins (sequences as provided by vendor)
extra <- tibble::tribble(
  ~DNA_ID,  ~Description,        ~Barcode,
  "C1374",  "CD159a/NKG2A",      "CAACTCCTGGGACTT",
  "C0054",  "CD34",              "GCAGAAATCTCCCTT",
  "C0061",  "CD117/c-kit",       "AGACTAATAGCTGAC"
)

all_abs <- dplyr::bind_rows(base, extra)

# Build feature_ref in the EXACT required order
feature_ref <- all_abs %>%
  dplyr::transmute(
    id           = make.unique(gsub("[^A-Za-z0-9_\\-]+", "_", DNA_ID)),
    name         = gsub("\\s+", "_", Description),  # name must not contain spaces
    read         = "R2",
    pattern      = "5PNNNNNNNNNN(BC)",
    sequence     = toupper(Barcode),                # forward (no reverse complement)
    feature_type = "Antibody Capture"
  ) %>%
  dplyr::distinct(sequence, .keep_all = TRUE) %>%
  dplyr::distinct(id, .keep_all = TRUE)

out <- "/net/beegfs/users/P089417/MDS_CRfiles/feature_ref.csv"
readr::write_csv(feature_ref, out)
message("Wrote feature_ref.csv with ", nrow(feature_ref), " rows at: ", out)

##Now this feature_ref.csv can directly be used as an input for CellRanger Multi
