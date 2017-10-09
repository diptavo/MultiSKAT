Open_SSD <- function(File.SSD, File.Info){
	return(SKAT::Open_SSD(File.SSD, File.Info))}

Close_SSD <- function(){
	return(SKAT::Close_SSD())}

Generate_SSD_SetID <- function(File.Bed, File.Bim, File.Fam, File.SetID,File.SSD, File.Info, Is.FlipGenotype=TRUE){
	return(SKAT::Generate_SSD_SetID(File.Bed, File.Bim, File.Fam, File.SetID,File.SSD, File.Info, Is.FlipGenotype=TRUE))}

Get_Genotypes_SSD <- function(SSD_INFO, Set_Index, is_ID = FALSE){
	return(SKAT::Get_Genotypes_SSD(SSD_INFO, Set_Index, is_ID = FALSE))}

Read_Plink_FAM <- function(Filename, Is.binary=TRUE, flag1=0){
	return(SKAT::Read_Plink_FAM(Filename, Is.binary=TRUE, flag1=0))}

Read_Plink_FAM_Cov <- function(Filename, File_Cov, Is.binary=TRUE, flag1=0, cov_header=TRUE){
	return(SKAT::Read_Plink_FAM_Cov(Filename, File_Cov, Is.binary=TRUE, flag1=0, cov_header=TRUE))}

Read_SNP_WeightFile <- function(FileName){
	return(SKAT::Read_SNP_WeightFile(FileName))}



