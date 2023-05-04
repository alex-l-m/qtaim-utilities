library(tidyverse)
library(glue)

indir <- commandArgs(trailingOnly = TRUE)[1]

paths <- tibble(path = Sys.glob(glue("{indir}/*_critpoints.csv"))) |>
    mutate(mol_id = str_match(basename(path), "(.*)_critpoints.csv")[,2])

mol_ids <- paths |>
  select(mol_id)

cp <- paths |>
    group_by(mol_id) |>
    reframe(read_csv(path, col_types = cols(
    critical_point = col_character(),
    property = col_character(),
    value = col_character())))

atom_tbl <- cp |>
    filter(property == "type_string" & value == "NACP") |>
    select(mol_id, critical_point) |>
    left_join(cp, by = c("mol_id", "critical_point")) |>
    filter(property == "atom_1") |>
    rename(atom_id = value) |>
    mutate(symbol = str_extract(atom_id, "[A-Z][a-z]*")) |>
    select(mol_id, critical_point, atom_id, symbol) |>
    left_join(cp, by = c("mol_id", "critical_point")) |>
    filter(str_detect(property, "Coords_[xyz]")) |>
    mutate(value = as.double(value)) |>
    mutate(coordinate = str_match(property, "Coords_([xyz])")[,2]) |>
    pivot_wider(id_cols = c(mol_id, critical_point, atom_id, symbol),
                names_from = coordinate, values_from = value) |>
    mutate(formal_charge = 0)

bond_tbl <- cp |>
    filter(property == "type_string" & value == "BCP") |>
    select(mol_id, critical_point) |>
    left_join(cp, by = c("mol_id", "critical_point")) |>
    filter(property == "atom_1" | property == "atom_2") |>
    pivot_wider(id_cols = c(mol_id, critical_point),
                names_from = property, values_from = value) |>
    rename(start_atom = atom_1, end_atom = atom_2) |>
    mutate(bond_id = sprintf("BCP%s", critical_point),
           bond_type = "UNSPECIFIED") |>
    select(mol_id, critical_point, bond_id, start_atom, end_atom, bond_type)

write_csv(mol_ids, "critpoint_mol_tbl.csv.gz")
write_csv(atom_tbl, "critpoint_one_tbl.csv.gz")
write_csv(bond_tbl, "critpoint_two_tbl.csv.gz")
