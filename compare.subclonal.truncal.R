rm(list = ls())
library(dplyr)
library(tidyr)
library(magrittr)

# Get MutationTime output
mt_output <- data.table::fread("data-mutation-timing-clonality.csv")

mt_output %>%
  dplyr::select(Patient.ID, Gene.mut, Truncal.branch, Clonality) %>%
  distinct() -> t1

# Make sure each mutation in a sector is always either Truncal or Branch but not both
stopifnot((t1 %>% 
  dplyr::select(Patient.ID, Gene.mut, Truncal.branch) %>%
  distinct() %>%
  group_by(Patient.ID, Gene.mut) %>%
  summarize(count = n(),.groups = "keep") %>%
  filter(count > 1) %>% nrow()) == 0)

# Create t2 with only the truncal mutations and with separte columns 
# for each possible value of Clonality
t1 %>%                     
  dplyr::filter(Truncal.branch == "Truncal") %>%
  distinct() %>%
  mutate(count = 1) %>%
  tidyr::pivot_wider(names_from = Clonality, values_from = count) -> t2

# Confirm that each mutation occurs only once for each patient
stopifnot((t2 %>% select(Patient.ID, Gene.mut) %>% distinct() %>% nrow()) ==
            nrow(t2))

# Put total truncal mutation counts for each patient
# in new table t2.total
t2 %>%
  group_by(Patient.ID) %>%
  summarize(total.count = n()) -> t2.total

# Put count of truncal mutations that are at least sometimes
# subclonal for each patient in  t2.subclonal
t2 %>%
  filter(!is.na(subclonal)) %>%
  group_by(Patient.ID) %>%
  summarize(subclonal.count = n()) -> t2.subclonal
           
left_join(t2.total, t2.subclonal) %>%
  mutate(subclonal.count = 
             ifelse(is.na(subclonal.count), 0, subclonal.count)) %>%
  mutate(patient_id = Patient.ID, 
         Patient.ID = NULL, 
         mt.clonal.count = total.count - subclonal.count) -> 
  mutation.timer.summary

t2 %>%
  filter(is.na(subclonal) & is.na(NA)) %>%
  select(Patient.ID, Gene.mut) %>%
  mutate(patient_id = Patient.ID, 
         mutation = Gene.mut, 
         Patient.ID = NULL, 
         Gene.mut = NULL) %>%
  group_by(patient_id) %>%
  arrange(mutation, .by_group = TRUE) -> 
  mutation.timer.clonal.truncal.mutations

# Get PyClone results
pc <- data.table::fread("pyclone/clonal.mut.pyclone.out.csv")

# Confirm that for each patients ID all mutations are
# have rows for the same number of sectors.
stopifnot((pc %>%
  group_by(patient_id, mutation_id_clonal) %>%
  summarise(sector.count = n()) %>%
  select(patient_id, sector.count) %>%
  distinct() %>% 
  group_by(patient_id) %>%
  summarise(count.sectors = n()) %>%
  filter(count.sectors > 1) %>%
  nrow()) == 0)

pc %>% 
  select(patient_id, mutation_id_clonal) %>%
  distinct() %>%
  mutate(mutation = mutation_id_clonal, mutation_id_clonal = NULL)-> 
  pc.clonal.truncal.mutations

pc.clonal.truncal.mutations %>%
  group_by(patient_id) %>%
  summarise(pc.clonal.count = n()) -> 
  pc.summary

full_join(mutation.timer.summary, pc.summary) %>%
  mutate(pc.subclonal.count = total.count - pc.clonal.count) ->
  compare

data.table::fwrite(compare, "mt.vs.pc.subclonal.truncal.csv")
