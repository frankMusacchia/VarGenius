
########################################PART 1
#The first part of this script must be run from Pico where all the files are present
#to create the table from which the heatmap will be created
# optionally you can donwload all the bed files and do it from your Desktop pc
#/cineca/prod/applications/r/3.2.2/gnu--4.8.3/bin/R

#COMMON VARIABLES
suffix = "baserecal_precalread_exonsnoncov"


#####VARIABLES SEPCIFIC FOR ANALYSIS
analtype = "trio"
gen_cov_suffix = "baserecal_precalread_Agilent_SureSelect_CRE_v2_UCSC_coding_exons.bed"
#panel = paste(analysisname,"final_gene_list",sep=".")
#CREv1
analyses <- c('UD_GE001','UD_GE002','UD_GE003','UD_GE004','UD_GE005','UD_GE006','UD_GE007','UD_GE008','UD_GE009','UD_GE010','UD_GE011','UD_GE012','UD_GE013','UD_GE015','UD_GE016','UD_GE017','UD_GE018','UD_GE019','UD_GE020','UD_GE021','UD_GE022','UD_GE023','UD_GE024','UD_GE025','UD_GE026','UD_GE027','UD_GE029','UD_GE030','UD_GE031','UD_GE032','UD_GE033','UD_GE034','UD_GE035','UD_GE036','UD_GE046','UD_MO002','UD_MO006','UD_MO008','UD_MO009','UD_MO015','UD_MO016','UD_MO017','UD_MO018','UD_MO019','UD_MO020','UD_MO021','UD_MO024','UD_MO027','UD_MO029','UD_MO030','UD_MO032','UD_MO033','UD_MO035','UD_MO036','UD_MO037','UD_MO038','UD_MO039','UD_MO040','UD_MO041','UD_MO043','UD_MO045','UD_MO046','UD_MO049','UD_MO050','UD_MO051','UD_MO055','UD_MO058','UD_MO062','UD_MO063','UD_MO068','UD_MO2016','UD_NA001','UD_NA002','UD_NA003','UD_NA004','UD_NA005','UD_NA008','UD_NA010','UD_NA012','UD_NA013','UD_NA014','UD_NA017','UD_NA018','UD_NA021','UD_NA022','UD_NA025','UD_NA026','UD_NA027','UD_NA028','UD_NA029','UD_NA030','UD_NA033','UD_NA036','UD_NA039','UD_NA040','UD_NA043','UD_NA045','UD_NA047','UD_NA050','UD_NA053','UD_NA054','UD_NA056','UD_NA057','UD_NA058','UD_NA059','UD_RM004','UD_RM005','UD_RM008','UD_RM013','UD_RM015','UD_RM1001','UD_RM1002','UD_RM1004','UD_RM1005','UD_RM1010','UD_RM1012','UD_RM1013','UD_RM1014','UD_RM1017','UD_RM1018','UD_RM1020','UD_RM1021','UD_X4111','UD_ZO002','UD_ZOL001','UD_ZOL004','UD_ZOL005','UD_ZOL008','UD_ZOL009','UD_ZOL010','UD_ZOL011','UD_ZOL012','UD_ZOL013','UD_ZOL014','UD_ZOL015')
#CREv2
analyses <- c('UD_BES001','UD_BES003','UD_BES004','UD_BES005','UD_BES006','UD_BES008','UD_BES010','UD_BES013','UD_BES014','UD_BES015','UD_BES016','UD_BES017','UD_BES018','UD_BES020','UD_BES021','UD_BES022','UD_BES023','UD_BES024','UD_BES026','UD_BES027','UD_BES028','UD_BES030','UD_GE039','UD_GE040','UD_GE041','UD_GE042','UD_GE043','UD_GE044','UD_GE045','UD_GE047','UD_GE050','UD_GE051','UD_GE052','UD_GE053','UD_GE054','UD_GE057','UD_GE058','UD_GE059','UD_GE060','UD_GE061','UD_GE062','UD_GE063','UD_GE064','UD_GE070','UD_GE071','UD_ME007','UD_ME010','UD_ME012','UD_MI002','UD_MO022','UD_MO031','UD_MO044','UD_MO048','UD_MO052','UD_MO053','UD_MO057','UD_MO059','UD_MO060','UD_MO061','UD_MO064','UD_MO065','UD_MO066','UD_MO069','UD_MO070','UD_MO072','UD_MO073','UD_MO074','UD_MO075','UD_MO076','UD_MO077','UD_MO078','UD_MO079','UD_MO080','UD_MO081','UD_MO082','UD_MO083','UD_MO084','UD_MO085','UD_MO086','UD_MO087','UD_MO089','UD_MO090','UD_MO091','UD_MO092','UD_MO093','UD_MO094','UD_MO095','UD_MO096','UD_MO097','UD_MO098','UD_MO099','UD_MO101','UD_MO103','UD_MO104','UD_MO105','UD_MO106','UD_MO107','UD_MO108','UD_MO109','UD_MO110','UD_MO111','UD_MO115','UD_MO117','UD_MO119','UD_MO122','UD_NA035','UD_NA041','UD_NA046','UD_NA048','UD_NA055','UD_NA062','UD_NA063','UD_NA064','UD_NA065','UD_NA066','UD_NA068','UD_NA069','UD_NA070','UD_NA071','UD_NA072','UD_NA074','UD_NA075','UD_NA076','UD_NA077','UD_NA078','UD_NA079','UD_NA081','UD_NA082','UD_NA083','UD_NA085','UD_NA087','UD_NA088','UD_NA089','UD_NA090','UD_NA091','UD_NA092','UD_NA095','UD_NA097','UD_NA098','UD_NA099','UD_NA100','UD_NA101','UD_NA102','UD_NA103','UD_NA104','UD_NA105','UD_NA106','UD_NA107','UD_NA108','UD_NA110','UD_NA111','UD_NA112','UD_OA003','UD_OA006','UD_OA007','UD_OA012','UD_RM1006','UD_RM1007','UD_RM1015','UD_RM1016','UD_RM1023','UD_RM1024','UD_RM1027','UD_RM1029','UD_RM1030','UD_RM1032','UD_RM1033','UD_RM1036','UD_RM1037','UD_RM1041','UD_RM1043','UD_RM1044','UD_RM1047','UD_ZOL018','UD_ZOL024','UD_ZOL025','UD_ZOL029')


#Mendeliomi Banfi
analtype = "single"
gen_cov_suffix="baserecal_precalread_Agilent_ClearSeq_inherited_diseases_exons_UCSC.bed"
panel = "Retinopathy"	
candpanel = "RetinopathyCAND"
analyses <- c('SB_ID_A174','SB_ID_A237','SB_ID_A354','SB_ID_A356','SB_ID_A358','SB_ID_A359','SB_ID_A362','SB_ID_A364','SB_ID_A375','SB_ID_A376','SB_ID_A377','SB_ID_A387','SB_ID_A391','SB_ID_A393_P1','SB_ID_A393_P2','SB_ID_A404_P1','SB_ID_A404_P2','SB_ID_A408','SB_ID_A414','SB_ID_A415','SB_ID_A416','SB_ID_A420','SB_ID_A423','SB_ID_A425','SB_ID_A427','SB_ID_A433','SB_ID_A434','SB_ID_A436','SB_ID_A437','SB_ID_A439_M','SB_ID_A439_P','SB_ID_A440','SB_ID_A443','SB_ID_A444','SB_ID_A445','SB_ID_A446','SB_ID_A447','SB_ID_A448','SB_ID_A449','SB_ID_A451','SB_ID_A452_rip','SB_ID_A453','SB_ID_A454','SB_ID_A456','SB_ID_A458_P','SB_ID_A459_F','SB_ID_A460','SB_ID_A474','SB_ID_A482_P','SB_ID_A483_M','SB_ID_A486','SB_ID_A498','SB_ID_A499_P','SB_ID_A501_F','SB_ID_A502','SB_ID_A509','SB_ID_A520','SB_ID_A523','SB_ID_A550','SB_ID_A551','SB_ID_A552','SB_ID_A553','SB_ID_A560','SB_ID_A561','SB_ID_A562','SB_ID_A563','SB_ID_A564','SB_ID_A566_P','SB_ID_A567_F','SB_ID_A572_F','SB_ID_A572_P','SB_ID_A590','SB_ID_A591','SB_ID_A592','SB_ID_A593','SB_ID_A612','SB_ID_A614_M','SB_ID_A614_P','SB_ID_A615_P1','SB_ID_A615_P2','SB_ID_A617','SB_ID_A621','SB_ID_A622','SB_ID_A623','SB_ID_A624','SB_ID_A625','SB_ID_A626','SB_ID_A646','SB_ID_A647','SB_ID_A654','SB_ID_A661','SB_ID_A665','SB_ID_A669','SB_ID_A688','SB_ID_A725','SB_ID_A737','SB_ID_A739','SB_ID_A740','SB_ID_A743','SB_ID_A744','SB_ID_A747','SB_ID_A748','SB_ID_A752','SB_ID_A753','SB_ID_A757','SB_ID_A785','SB_ID_A788','SB_ID_A791','SB_ID_A792','SB_ID_A793','SB_ID_A794','SB_ID_A795','SB_ID_A846','SB_ID_A855','SB_ID_A860','SB_ID_A864','SB_ID_A902','SB_ID_A914','SB_ID_A922','SB_ID_A928','SB_ID_A934','SB_ID_A954','SB_ID_A955','SB_ID_A956','SB_ID_A958','SB_ID_A964','SB_ID_A966','SB_ID_A976','SB_ID_ARRP133','SB_ID_ARRP158','SB_ID_ARRP162','SB_ID_ARRP165','SB_ID_ARRP182','SB_ID_USH197')

#CCP Banfi mancano 'SB_CCP_A862',SB_CCP_1998-17,'SB_CCP_A945',
gen_cov_suffix="baserecal_precalread_Agilent_SureSelect_CCP_v1_exons_UCSC.bed"
analyses <- c('SB_CCP_1998-17','SB_CCP_A211','SB_CCP_A286','SB_CCP_A301','SB_CCP_A555','SB_CCP_A559','SB_CCP_A587','SB_CCP_A603','SB_CCP_A639','SB_CCP_A666','SB_CCP_A681','SB_CCP_A704','SB_CCP_A852','SB_CCP_A918','SB_CCP_A929','SB_CCP_A937','SB_CCP_A951','SB_CCP_A959','SB_CCP_X5185','SB_CCP_X5186','SB_CCP_X5187','SB_CCP_X5188','SB_CCP_X5410_P','SB_CCP_X5410_F','SB_CCP_X5410_M','SB_CCP_X5452','SB_CCP_X5457','SB_CCP_X5672','VN_CCP_X4925','VN_CCP_X5023','VN_CCP_X5024','VN_CCP_X5094','VN_CCP_X5095','VN_CCP_X5096','VN_CCP_X5097','VN_CCP_X5098','VN_CCP_X5099','VN_CCP_X5100','VN_CCP_X5107','VN_CCP_X5157','VN_CCP_X5159','VN_CCP_X5160','VN_CCP_X5161','VN_CCP_X5168','VN_CCP_X5169','VN_CCP_X5170','VN_CCP_X5318')

#CREv1 missing: ,'SB_ClinExome_A292'
gen_cov_suffix="baserecal_precalread_Agilent_SureSelect_CRE_v1_exons_UCSC.bed"

analyses <- c('SB_A217','SB_A234','SB_A252','SB_A254','SB_A268','SB_A290','SB_A313','SB_A314','SB_A325','SB_A327','SB_A389_P','SB_A390','SB_A392','SB_ARRP129','SB_ClinExome_A314','SB_ClinExome_A337','SB_ClinExome_A348','SB_ClinExome_A368','SB_ClinExome_A374','SB_ClinExome_A395','SB_WE_S4100_F','SB_WE_S4100_M','SB_WE_S4100_P1','SB_WE_S4100_P2')

#Nigro 
#ID
gen_cov_suffix = "baserecal_precalread_Agilent_ClearSeq_inherited_diseases_exons_UCSC.bed"
#trios
analtype = "trio"
panel = paste(analysisname,"final_gene_list",sep=".")
analyses <- c('ID_L363','ID_NA037','ID_NA041','ID_RM1007','ID_X2596','ID_X3618','ID_X4007','ID_X4071','ID_X4118','ID_X4130','ID_X4169','ID_X4218','ID_X4221','ID_X4288','ID_X4315','VN_ID_NA048','VN_ID_NA055','VN_ID_X3725','VN_ID_X3858','VN_ID_X4350','VN_ID_X4368','VN_ID_X4375','VN_ID_X4378','VN_ID_X4387','VN_ID_X4390','VN_ID_X4393','VN_ID_X4398','VN_ID_X4405','VN_ID_X4417','VN_ID_X4423','VN_ID_X4429','VN_ID_X4432_II','VN_ID_X4451','VN_ID_X4454','VN_ID_X4460','VN_ID_X4464','VN_ID_X4467','VN_ID_X4470','VN_ID_X4473','VN_ID_X4479','VN_ID_X4484','VN_ID_X4490','VN_ID_X4492','VN_ID_X4501','VN_ID_X4504','VN_ID_X4516','VN_ID_X4537','VN_ID_X4552','VN_ID_X4566','VN_ID_X4570','VN_ID_X4585','VN_ID_X4613','VN_ID_X4616','VN_ID_X4619','VN_ID_X4623','VN_ID_X4627','VN_ID_X4630','VN_ID_X4648','VN_ID_X4658','VN_ID_X4666','VN_ID_X4669_RIP','VN_ID_X4676','VN_ID_X4679','VN_ID_X4685','VN_ID_X4695','VN_ID_X4698','VN_ID_X4703','VN_ID_X4711','VN_ID_X4714','VN_ID_X4720','VN_ID_X4723','VN_ID_X4738','VN_ID_X4768','VN_ID_X4771','VN_ID_X4774','VN_ID_X4780','VN_ID_X4878','VN_ID_X4889','VN_ID_X4892','VN_ID_X4919','VN_ID_X4927','VN_ID_X4945','VN_ID_X4953','VN_ID_X4966','VN_ID_X4970','VN_ID_X4987','VN_ID_X5010','VN_ID_X5013','VN_ID_X5016','VN_ID_X5019','VN_ID_X5035','VN_ID_X5038','VN_ID_X5061','VN_ID_X5070','VN_ID_X5073','VN_ID_X5076','VN_ID_X5171','VN_ID_X5174','VN_ID_X5197')

#singoli
analtype = "single"
#analyses <- c('ID_K91_P','ID_NA037_P','ID_X4128_P','ID_X4179_P','ID_X4219_P','VN_ID_X4365_P','VN_ID_X4556_P','VN_ID_X4578_P','VN_ID_X4594_P','VN_ID_X4717_II_P','VN_ID_X4726_P','VN_ID_X4735_II_P','VN_ID_X4913_P','ID_NA044_P','VN_ID_X4118_P','VN_ID_X4706_P','VN_ID_X4913_P','ID_K90','ID_L369','ID-NA016','ID_X1496','ID_X2122','ID_X3065','ID_X3359','ID_X3791','ID_X3792','ID_X3803','ID_X3830','ID_X3848','ID_X3849','ID_X3866','ID_X3894','ID-X3936','ID-X3937','ID-X3939','ID-X3941','ID-X3942','ID-X3943','ID-X3986','ID_X3987','ID_X3998','ID_X3999','ID-X4008','ID_X4038','ID_X4155','ID_X4156','ID_X4183','ID_X4217','ID_X4253','ID_X4287','ID_X4296','ID_X4318','VN_ID_DMD10368','VN_ID_DMD105','VN_ID_DMD119','VN_ID_DMD130','VN_ID_DMD211','VN_ID_DMD303','VN_ID_DMD3201','VN_ID_DMD3261','VN_ID_DMD419','VN_ID_DMD463','VN_ID_DMD4753','VN_ID_DMD4813','VN_ID_DMD5160','VN_ID_DMD5514','VN_ID_DMD813','VN_ID_HF74','VN_ID_K120','VN_ID_NF759','VN_ID_NF760','VN_ID_NF762','VN_ID_NF768','VN_ID_NF769','VN_ID_NF770','VN_ID_NF771','VN_ID_X3468','VN_ID_X4123','VN_ID_X4157','VN_ID_X4213','VN_ID_X4227_II','VN_ID_X4276','VN_ID_X4284','VN_ID_X4327','VN_ID_X4355','VN_ID_X4367','VN_ID_X4373','VN_ID_X4374','VN_ID_X4381','VN_ID_X4385','VN_ID_X4396','VN_ID_X4397','VN_ID_X4401','VN_ID_X4402','VN_ID_X4414','VN_ID_X4415','VN_ID_X4416','VN_ID_X4418','VN_ID_X4428','VN_ID_X4435','VN_ID_X4436','VN_ID_X4440','VN_ID_X4444','VN_ID_X4445','VN_ID_X4446','VN_ID_X4447','VN_ID_X4487','VN_ID_X4488','VN_ID_X4489','VN_ID_X4495','VN_ID_X4496','VN_ID_X4498','VN_ID_X4499','VN_ID_X4500','VN_ID_X4515','VN_ID_X4519','VN_ID_X4542','VN_ID_X4543','VN_ID_X4544','VN_ID_X4545','VN_ID_X4546','VN_ID_X4547','VN_ID_X4549','VN_ID_X4555','VN_ID_X4558','VN_ID_X4562','VN_ID_X4563','VN_ID_X4564','VN_ID_X4565','VN_ID_X4569','VN_ID_X4584','VN_ID_X4593','VN_ID_X4596','VN_ID_X4597','VN_ID_X4598','VN_ID_X4606','VN_ID_X4607','VN_ID_X4610','VN_ID_X4616','VN_ID_X4617','VN_ID_X4619','VN_ID_X4620','VN_ID_X4622','VN_ID_X4626','VN_ID_X4633','VN_ID_X4636','VN_ID_X4639','VN_ID_X4646','VN_ID_X4651','VN_ID_X4655','VN_ID_X4663','VN_ID_X4664','VN_ID_X4665','VN_ID_X4672','VN_ID_X4674','VN_ID_X4675','VN_ID_X4701','VN_ID_X4704','VN_ID_X4737','VN_ID_X4741','VN_ID_X4750','VN_ID_X4751','VN_ID_X4757','VN_ID_X4777','VN_ID_X4792','VN_ID_X4795','VN_ID_X4895','VN_ID_X4896','VN_ID_X4897','VN_ID_X4912','VN_ID_X4918','VN_ID_X4922','VN_ID_X4923','VN_ID_X4924','VN_ID_X4925','VN_ID_X4926','VN_ID_X4936','VN_ID_X4940','VN_ID_X4943','VN_ID_X4944','VN_ID_X4948','VN_ID_X4949','VN_ID_X4950','VN_ID_X4951','VN_ID_X4952','VN_ID_X4954','VN_ID_X4955','VN_ID_X4956','VN_ID_X4960','VN_ID_X4961','VN_ID_X4962','VN_ID_X4969','VN_ID_X4983','VN_ID_X4984','VN_ID_X4985','VN_ID_X4993','VN_ID_X4994','VN_ID_X4995','VN_ID_X4996','VN_ID_X4997','VN_ID_X4998','VN_ID_X5001','VN_ID_X5002','VN_ID_X5022','VN_ID_X5023','VN_ID_X5024','VN_ID_X5025','VN_ID_X5041','VN_ID_X5044','VN_ID_X5045','VN_ID_X5048','VN_ID_X5049','VN_ID_X5050','VN_ID_X5051','VN_ID_X5052','VN_ID_X5053','VN_ID_X5054','VN_ID_X5055','VN_ID_X5060','VN_ID_X5064','VN_ID_X5067','VN_ID_X5147','VN_ID_X5148','VN_ID_X5149','VN_ID_X5156')
analyses <- c('VN_ID_NF759','VN_ID_NF760','VN_ID_NF769','VN_ID_NF770','VN_ID_NF771','VN_ID_X3468','VN_ID_X3868_F','VN_ID_X3868_M','VN_ID_X3868_P1','VN_ID_X3868_P2','VN_ID_X3868_P3','VN_ID_X4617','VN_ID_X4630_F','VN_ID_X4630_P1','VN_ID_X4630_P2','VN_ID_X4717_II_P','VN_ID_X4723_F','VN_ID_X4723_M','VN_ID_X4723_P','VN_ID_X4735_II_P','VN_ID_X4737','VN_ID_X4738_M','VN_ID_X4738_P1','VN_ID_X4738_P2','VN_ID_X4741','VN_ID_X4750','VN_ID_X4777','VN_ID_X4889_F','VN_ID_X4889_M','VN_ID_X4889_P','VN_ID_X4897','VN_ID_X4913_P1','VN_ID_X4913_P2','VN_ID_X4918','VN_ID_X4919_F','VN_ID_X4919_M','VN_ID_X4919_P','VN_ID_X4923','VN_ID_X4924','VN_ID_X4925','VN_ID_X4926','VN_ID_X4927_F','VN_ID_X4927_M','VN_ID_X4927_P','VN_ID_X4936','VN_ID_X4940','VN_ID_X4948','VN_ID_X4949','VN_ID_X4952','VN_ID_X4955','VN_ID_X4966_F','VN_ID_X4966_M','VN_ID_X4966_P','VN_ID_X4969','VN_ID_X4985','VN_ID_X4994','VN_ID_X4995','VN_ID_X4996','VN_ID_X4998','VN_ID_X5002_P','VN_ID_X5013_F','VN_ID_X5013_M','VN_ID_X5013_P','VN_ID_X5016_F','VN_ID_X5016_M','VN_ID_X5016_P','VN_ID_X5022','VN_ID_X5023','VN_ID_X5024','VN_ID_X5025','VN_ID_X5038_F','VN_ID_X5038_M','VN_ID_X5038_P','VN_ID_X5045','VN_ID_X5048','VN_ID_X5051','VN_ID_X5052','VN_ID_X5055','VN_ID_X5060','VN_ID_X5070_F','VN_ID_X5070_M','VN_ID_X5070_P','VN_ID_X5073_F','VN_ID_X5073_M','VN_ID_X5073_P','VN_ID_X5147','VN_ID_X5149','VN_ID_X5197_F','VN_ID_X5197_M','VN_ID_X5197_P')

#CREv1 - NON SONO IN VARGENIUSANALYSES
gen_cov_suffix = "baserecal_precalread_Agilent_SureSelect_CRE_v1_UCSC_coding_exons.bed"
#TRIOS
analtype = "trio"
panel = paste(analysisname,"final_gene_list",sep=".")
analyses <- c('VN_CCP_NF817_TRIO','VN_CCP_X4378_TRIO','VN_CCP_X4558_TRIO','VN_CCP_X4789_TRIO','VN_CCP_X4796_TRIO','VN_CCP_X4799_TRIO','VN_CCP_X5004_TRIO','VN_CCP_X5007_TRIO','VN_CCP_X5102_TRIO','VN_CCP_X5110_TRIO','VN_CCP_X5129_TRIO','VN_CCP_X5136_TRIO','VN_CCP_X5139_TRIO','VN_CCP_X5153_TRIO','VN_CCP_X5178_TRIO','VN_CCP_X5191_TRIO','VN_CCP_X5250_TRIO','VN_CCP_X5274_TRIO','VN_CCP_X5283_TRIO','VN_CCP_X5286_TRIO','VN_CCP_X5289_TRIO','VN_CCP_X5310_TRIO','VN_CCP_X5320_TRIO','VN_CCP_X5383_TRIO','VN_CCP_X5433_TRIO','VN_CCP_X5444_TRIO','VN_CCP_X5559_TRIO','VN_CCP_X5650_TRIO')
#Singoli
analtype = "single"
gen_cov_suffix="baserecal_precalread_Agilent_SureSelect_CCP_v1_exons_UCSC.bed"
#analyses <- c('VN_CCP_NF817_P','VN_CCP_X4378_P','VN_CCP_X4558_P','VN_CCP_X4789_P','VN_CCP_X4796_P','VN_CCP_X4799_P','VN_CCP_X5004_P','VN_CCP_X5007_P','VN_CCP_X5102_P','VN_CCP_X5110_P','VN_CCP_X5129_P','VN_CCP_X5136_P','VN_CCP_X5139_P','VN_CCP_X5153_P','VN_CCP_X5178_P','VN_CCP_X5191_P','VN_CCP_X5250_P','VN_CCP_X5274_P','VN_CCP_X5283_P','VN_CCP_X5286_P','VN_CCP_X5289_P','VN_CCP_X5310_P','VN_CCP_X5320_P','VN_CCP_X5383_P','VN_CCP_X5433_P','VN_CCP_X5444_P','VN_CCP_X5559_P','VN_CCP_X5650_P','VN_CCP_NF817_M','VN_CCP_X4378_M','VN_CCP_X4558_M','VN_CCP_X4789_M','VN_CCP_X4796_M','VN_CCP_X4799_M','VN_CCP_X5004_M','VN_CCP_X5007_M','VN_CCP_X5102_M','VN_CCP_X5110_M','VN_CCP_X5129_M','VN_CCP_X5136_M','VN_CCP_X5139_M','VN_CCP_X5153_M','VN_CCP_X5178_M','VN_CCP_X5191_M','VN_CCP_X5250_M','VN_CCP_X5274_M','VN_CCP_X5283_M','VN_CCP_X5286_M','VN_CCP_X5289_M','VN_CCP_X5310_M','VN_CCP_X5320_M','VN_CCP_X5383_M','VN_CCP_X5433_M','VN_CCP_X5444_M','VN_CCP_X5559_M','VN_CCP_X5650_M','VN_CCP_NF817_F','VN_CCP_X4378_F','VN_CCP_X4558_F','VN_CCP_X4789_F','VN_CCP_X4796_F','VN_CCP_X4799_F','VN_CCP_X5004_F','VN_CCP_X5007_F','VN_CCP_X5102_F','VN_CCP_X5110_F','VN_CCP_X5129_F','VN_CCP_X5136_F','VN_CCP_X5139_F','VN_CCP_X5153_F','VN_CCP_X5178_F','VN_CCP_X5191_F','VN_CCP_X5250_F','VN_CCP_X5274_F','VN_CCP_X5283_F','VN_CCP_X5286_F','VN_CCP_X5289_F','VN_CCP_X5310_F','VN_CCP_X5320_F','VN_CCP_X5383_F','VN_CCP_X5433_F','VN_CCP_X5444_F','VN_CCP_X5559_F','VN_CCP_X5650_F','VN_CCP_X5003_P','VN_CCP_X5108_P','VN_CCP_X5134_P','VN_CCP_X5216_P','VN_ID_X5002_P','VN_CCP_X5113_P','VN_CCP_X5163_P','VN_CCP_X5392_P','VN_CCP_A909','VN_CCP_A930','VN_CCP_A965','VN_CCP_NF750','VN_CCP_NF769','VN_CCP_NF772','VN_CCP_NF774','VN_CCP_NF775','VN_CCP_NF779','VN_CCP_NF780','VN_CCP_NF781','VN_CCP_NF790','VN_CCP_NF792','VN_CCP_NF796','VN_CCP_NF802','VN_CCP_NF819','VN_CCP_SM207','VN_CCP_SM218','VN_CCP_X4682','VN_CCP_X4683','VN_CCP_X4684','VN_CCP_X4789','VN_CCP_X4794','VN_CCP_X4795','VN_CCP_X4965','VN_CCP_X4973','VN_CCP_X4998','VN_CCP_X5067','VN_CCP_X5086','VN_CCP_X5088','VN_CCP_X5089','VN_CCP_X5093','VN_CCP_X5101','VN_CCP_X5105','VN_CCP_X5106','VN_CCP_X5118','VN_CCP_X5119','VN_CCP_X5127','VN_CCP_X5132','VN_CCP_X5133','VN_CCP_X5142','VN_CCP_X5151','VN_CCP_X5152','VN_CCP_X5162','VN_CCP_X5167','VN_CCP_X5177','VN_CCP_X5181','VN_CCP_X5182','VN_CCP_X5183','VN_CCP_X5184','VN_CCP_X5404','VN_CCP_X5546','VN_CCP_X5747')
analyses <- c('VN_CCP_A909','VN_CCP_A930','VN_CCP_A965','VN_CCP_X2901_M','VN_CCP_X2901_P1','VN_CCP_X2901_P2','VN_CCP_X4378_F','VN_CCP_X4378_M','VN_CCP_X4378_P','VN_CCP_X4558_F','VN_CCP_X4558_M','VN_CCP_X4558_P','VN_CCP_X4682','VN_CCP_X4683','VN_CCP_X4684','VN_CCP_X4789_P','VN_CCP_X4794','VN_CCP_X4795','VN_CCP_X4796_F','VN_CCP_X4796_M','VN_CCP_X4796_P','VN_CCP_X4799_F','VN_CCP_X4799_M','VN_CCP_X4799_P','VN_CCP_X4965','VN_CCP_X4973','VN_CCP_X4998','VN_CCP_X5003_M','VN_CCP_X5004_F','VN_CCP_X5004_M','VN_CCP_X5004_P','VN_CCP_X5007_F','VN_CCP_X5007_M','VN_CCP_X5007_P','VN_CCP_X5067','VN_CCP_X5086','VN_CCP_X5088','VN_CCP_X5089','VN_CCP_X5093','VN_CCP_X5101','VN_CCP_X5102_F','VN_CCP_X5102_M','VN_CCP_X5102_P','VN_CCP_X5105','VN_CCP_X5106','VN_CCP_X5108_M','VN_CCP_X5108_P','VN_CCP_X5110_F','VN_CCP_X5110_M','VN_CCP_X5110_P','VN_CCP_X5113_F','VN_CCP_X5113_M','VN_CCP_X5113_P1','VN_CCP_X5113_P2','VN_CCP_X5118','VN_CCP_X5119','VN_CCP_X5127','VN_CCP_X5129_F','VN_CCP_X5129_M','VN_CCP_X5129_P','VN_CCP_X5132','VN_CCP_X5133','VN_CCP_X5134_M','VN_CCP_X5134_P','VN_CCP_X5136_P1','VN_CCP_X5136_P2','VN_CCP_X5136_P3','VN_CCP_X5139_F','VN_CCP_X5139_M','VN_CCP_X5139_P','VN_CCP_X5142','VN_CCP_X5151','VN_CCP_X5152','VN_CCP_X5153_F','VN_CCP_X5153_M','VN_CCP_X5153_P','VN_CCP_X5162','VN_CCP_X5163_F','VN_CCP_X5163_M','VN_CCP_X5163_P1','VN_CCP_X5163_P2','VN_CCP_X5167','VN_CCP_X5177','VN_CCP_X5178_F','VN_CCP_X5178_M','VN_CCP_X5178_P','VN_CCP_X5181','VN_CCP_X5182','VN_CCP_X5183','VN_CCP_X5184','VN_CCP_X5191_F','VN_CCP_X5191_M','VN_CCP_X5191_P','VN_CCP_X5216_F','VN_CCP_X5216_P','VN_CCP_X5250_F','VN_CCP_X5250_M','VN_CCP_X5250_P','VN_CCP_X5274_F','VN_CCP_X5274_M','VN_CCP_X5274_P','VN_CCP_X5283_F','VN_CCP_X5283_M','VN_CCP_X5283_P','VN_CCP_X5286_F','VN_CCP_X5286_M','VN_CCP_X5286_P','VN_CCP_X5289_F','VN_CCP_X5289_M','VN_CCP_X5289_P','VN_CCP_X5310_F','VN_CCP_X5310_M','VN_CCP_X5310_P','VN_CCP_X5320_F','VN_CCP_X5320_M','VN_CCP_X5320_P','VN_CCP_X5383_F','VN_CCP_X5383_M','VN_CCP_X5383_P','VN_CCP_X5392_F','VN_CCP_X5392_M','VN_CCP_X5392_P1','VN_CCP_X5392_P2','VN_CCP_X5404','VN_CCP_X5433_F','VN_CCP_X5433_M','VN_CCP_X5433_P','VN_CCP_X5444_F','VN_CCP_X5444_M','VN_CCP_X5444_P','VN_CCP_X5546','VN_CCP_X5559_F','VN_CCP_X5559_M','VN_CCP_X5559_P','VN_CCP_X5650_F','VN_CCP_X5650_M','VN_CCP_X5650_P','VN_CCP_X5747')
#####################################ANALYSIS 2

interval_2_gene <- data.frame(compid=character(),
                 Gene=character(),
                 stringsAsFactors=FALSE) 
                 
nc5g_all <- data.frame(Date=as.Date(character()),
                 File=character(), 
                 User=character(), 
                 stringsAsFactors=FALSE) 

#analysis_name="UD_BES002_P1"
analnum = 1
for (analysis in analyses){
		if (analtype == 'trio'){
		analysis_name = paste(analysis,"P",sep="_")
		}else{analysis_name = analysis}	
		nc5genes_t <- paste(analysis_name,suffix,sep="_")
		print (nc5genes_t)
		#Check if file is empty
		if (file.info(nc5genes_t)$size > 0){
			nc5genes <- read.delim(nc5genes_t,stringsAsFactors=F,head=F,quote="",sep="\t")

			
			#Get a composite id with chr_start_end of the exon
			nc5genes$compid <-  paste(nc5genes$V1,paste(nc5genes$V2,nc5genes$V3,sep="-"),sep=":")
			nc5genes_filt <- nc5genes[,c("compid","V7")]
			
			#Set column names
			colnames(nc5genes_filt) <- c("exon",analysis_name)
			
			
			if (nrow(nc5genes_filt) > 0){

				#Put it as numeric
				nc5genes_filt[[analysis_name]] <- as.numeric(nc5genes_filt[[analysis_name]])
				#Define the complete file if it is empty
				if (analnum == 1){
					nc5g_all <- nc5genes_filt
					interval_2_gene <- nc5genes[,c("V11","compid")]
					#colnames(genes_corresp) <- c("Chr","Start","End","Gene")
				}else{
					#Merge with the all file	
					nc5g_all <- merge(nc5g_all,nc5genes_filt,by="exon",all=T)
					nc5genes_2 <-nc5genes[!(nc5genes$compid %in% interval_2_gene$compid),c("V11","compid")]
					interval_2_gene <- rbind(interval_2_gene,nc5genes_2)
				}
				analnum = analnum + 1			
			}		
		}	
}



	suspect_table <- data.frame(Interval=character(),
					 samplename=character(), 
					 Value=character(), 
					 stringsAsFactors=FALSE) 
					 
	#Now count the occurrencies of percentages.  Equal percentages are most probably regions which cannot be
	#seen by the sequencer
	out_name = "suspect_intervals.txt"
	cat(paste("interv_uniq","sam_name_uniq","value",sep="\t"),file=out_name,append = FALSE)
	cat("\n",file=out_name,append = TRUE)

	interv_num = 0
	#Define the complete file if it is empty
	if (nrow(nc5g_all) > 0){
		nc5g_all[is.na(nc5g_all)] <- 1
		nc5g_all_t <- t(nc5g_all)
		#For each interval evaluate the proportion of coverage
		for ( interval in 2:ncol(nc5g_all_t)){
			#Count the number of time each percentage is present
			freq_table <-as.data.frame(ftable(nc5g_all_t[,interval]))
			#It counts also the interval name. Hence I remove this count
			freq_table_ok <- freq_table[1:nrow(freq_table)-1,]
			#Leave only frequencies which could be useful
			freq_table_ok_f <- freq_table_ok[as.vector(freq_table_ok$Var1) < 1 & as.vector(freq_table_ok$Freq) < 2,]
			
			#Only if there is only a single sample with a low coverage, go ahed
			if (nrow(freq_table_ok_f) == 1 ){
				#Get the value
				value <- as.numeric(as.vector(freq_table_ok_f$Var1)[1])
				#Get the name of the sample which has this "unique" low coverage
				samplename <- colnames(nc5g_all[which(nc5g_all[interval,]==value)])
				#Get the interval for which this sample has a SV
				interv_uniq <- nc5g_all$exon[interval]
	
				#Add  interval,samplename,value to a list of suspect SV
				suspect_table[nrow(suspect_table) + 1,] = list(interv_uniq,samplename,value)
			}
		}
	}	

	###Evaluation of suspect low-coverage intervals
	#The suspect non covered exons now should be compared with the other interval 
	#Those of the parents and also we need the average number of reads present in 
	#the other samples, beecause it can happen that other samples BAM cover the interval but
	#with very few reads	
	if ( nrow(suspect_table) > 0){
		
		#Get the compid (intervals) to be used
		uniq_int <- as.data.frame(unique(suspect_table$Interval))
		colnames(uniq_int) <- c("compid")
		geno_cov_all <- data.frame(Date=as.Date(character()),
						 File=character(), 
						 User=character(), 
						 stringsAsFactors=FALSE) 

		if (analtype == 'trio'){
			family_suffs = c("P","F","M")
		}else{family_suffs = c("P")}	
								 
		analnum = 1
		
		#Get the information for the suspect intervals from the genomeCoverageBed files for both probands and parents
		for (analysis_name in analyses){
			for (family_suff in family_suffs) {
				if (analtype == 'trio'){
					analysis = paste(analysis_name,family_suff,sep="_")
				}else{	analysis = analysis_name}	
		
				geno_cov_f <- paste(analysis,gen_cov_suffix,sep="_")
				#Check that the file exists and is non empty
				if (file.info(geno_cov_f)$size > 0 & file.exists(geno_cov_f)){
					geno_cov <- read.delim(geno_cov_f,stringsAsFactors=F,head=F,quote="",sep="\t")
					#Get a composite id with chr_start_end of the exon
					
					geno_cov$compid <- paste(geno_cov$V1,paste(geno_cov$V2,geno_cov$V3,sep="-"),sep=":")
					#Filter only the needed intervals
					geno_cov_needed <- merge(geno_cov,uniq_int,by="compid",all=F)
					#Select columns to investigate
					geno_cov_filt <- geno_cov_needed[,c("compid","V4")]
					colnames(geno_cov_filt) <- c("compid",analysis)
					#Define the complete file if it is empty
					if (analnum == 1){
						geno_cov_all <- geno_cov_filt
						
					}else{
						#Merge with the all file	
						geno_cov_all <- merge(geno_cov_all,geno_cov_filt,by="compid",all=T)
					}	
					analnum = analnum+1						
				}

			}
		}
		
		#Convert to numeric
		geno_cov_all_n<-as.matrix(sapply(geno_cov_all[, 2:ncol(geno_cov_all)], as.numeric))
		#re-Assign row names
		rownames(geno_cov_all_n)<-geno_cov_all$compid
		#Get the average for eachrow 
		geno_cov_all_n_st <- cbind(geno_cov_all_n,mean=rowMeans(geno_cov_all_n),compid=geno_cov_all$compid)

		#For each uniq interval generate a table with
		#samplename - interval value - prob_reads - fath_reads - moth_reads -avg_reads
		uniq_int<-as.matrix(uniq_int)
		
		##Extract data for TRIO
		if (analtype == 'trio'){
					
			final_table <- data.frame(
                 samplename=character(), 
                 compid=character(),
                 value=character(),
                  preads=character(),  
                  freads=character(),  
                  mreads=character(),  
                  avgreads=character(),  
                 stringsAsFactors=FALSE) 
			for (interv_num in 1:length(uniq_int) ){

				samplename <- suspect_table[suspect_table$Interval == uniq_int[interv_num],c("samplename")]
				list = unlist(strsplit(samplename,"_"))
				a<-list[1:length(list)-1]
				analysis <- paste(a,collapse='_')
				value <- suspect_table[suspect_table$Interval == uniq_int[interv_num],c("Value")]
				
				preads <- geno_cov_all[geno_cov_all$compid == uniq_int[interv_num],c(paste(analysis,"P",sep="_"))]	
				freads <- geno_cov_all[geno_cov_all$compid == uniq_int[interv_num],c(paste(analysis,"F",sep="_"))]	
				mreads <- geno_cov_all[geno_cov_all$compid == uniq_int[interv_num],c(paste(analysis,"M",sep="_"))]	
				
				if (length(preads) & length(freads) & length(mreads)){
					avgreads <- geno_cov_all_n_st[as.data.frame(geno_cov_all_n_st)$compid == uniq_int[interv_num],c("mean")]
					final_table[nrow(final_table) + 1,] = list(analysis,uniq_int[interv_num],value,preads,mreads,freads,avgreads)				
				}
			}		
		}else{	#Extract data for single samples
			final_table <- data.frame(
                 samplename=character(), 
                 compid=character(),
                 value=character(),
                  preads=character(), 
                  avgreads=character(),  
                 stringsAsFactors=FALSE) 
			for (interv_num in 1:length(uniq_int) ){

				samplename <- suspect_table[suspect_table$Interval == uniq_int[interv_num],c("samplename")]
				value <- suspect_table[suspect_table$Interval == uniq_int[interv_num],c("Value")]
				analysis <- samplename
				
				preads <- geno_cov_all[geno_cov_all$compid == uniq_int[interv_num],c(analysis)]	
				
				if (length(preads) ){
					avgreads <- geno_cov_all_n_st[as.data.frame(geno_cov_all_n_st)$compid == uniq_int[interv_num],c("mean")]
					final_table[nrow(final_table) + 1,] = list(analysis,uniq_int[interv_num],value,preads,avgreads)				
				}
			}	
		}	
				
		final_table <- final_table[order(final_table$samplename),]
		write.table(final_table,file="final_suspect_table.txt",sep="\t",row.names=F, quote=F)

	}

####SINGLE ANALYSIS FILTER and ADD GENE NAME

	#Count the number of time each percentage is present
	freq_table <- as.data.frame(ftable(final_table$samplename))
	freq_table_f <- freq_table[freq_table$Freq<100,]
	final_table_f <- merge(final_table,freq_table_f,by.x="samplename",by.y="Var1",all.x=F,all.y=F)
	
	#Add gene name using the interval gene correspondence
	final_table_fg <- merge(final_table_f,interval_2_gene,by="compid",all.x=T)
	write.table(final_table_fg,file="final_suspect_table_f.txt",sep="\t",row.names=F, quote=F)
	
	#If a panel of gene is present, add the column
	panelset="../gene_2_panel.txt"
	panel=""
	if ( panel != ""){
		tag_panel_genes("final_suspect_table_f.txt",panel,panel,candpanel,"-","V11",paste(gen_cov_suffix,"final_suspect_table_f_panel.txt",sep="_"))	
	}else{
		if ( panelset != "" ){
			g2p <- read.delim(panelset,stringsAsFactors=F,head=F,quote="",sep="\t")
			final_table_fgp <- merge(final_table_fg,g2p,by.x="V11",by.y="V1",all.x=T,all.y=F)
			write.table(final_table_fgp,file="final_suspect_table_gp.txt",sep="\t",row.names=F, quote=F)
		}
	}
	
#This function Given panel genes and candidate genes to belong to the panel, will add
#a column with the panel name and the field contains TRUE if the gene
#is one from the panel and CAND if belong to the candidates list
#Panel and candidate list must have column name called "Gene"
tag_panel_genes = function (table1,panel_name,panel_path,cand_path,empty_val,gene_field,out_name) {
	trueflag = "TRUE"
	candflag = "CAND"
	panel <- read.table(panel_path,stringsAsFactors=F,quote="",head=T,sep="\t")
	cand <- read.table(cand_path,stringsAsFactors=F,quote="",head=T,sep="\t")
	
	a<-read.delim(table1,stringsAsFactors=F,quote="",head=T,sep="\t")
	#Add the flag-column with the panel name
	a[[panel_name]] <-  ifelse(a[[gene_field]] %in% panel$Gene, trueflag, ifelse(a[[gene_field]] %in% cand$Gene, candflag, empty_val))
	#Print out a table
	write.table(a,file=out_name,sep="\t",quote=F,row.names=F)
}

############################################# ANALYSIS OF HETEROZYGOUS CNVs and DUPLICATIONS

#For this analysis I use the output from coverageBed reporting:
# The number of features in B (reads) that overlapped (by at least one base pair) the A interval.
# The number of bases in A that had non-zero coverage from features in B.
# The length of the entry in A.
# The fraction of bases in A that had non-zero coverage from features in B.

geno_cov_all <- data.frame(Date=as.Date(character()),
						 File=character(), 
						 User=character(), 
						 stringsAsFactors=FALSE) 

if (analtype == 'trio'){
	family_suffs = c("P","F","M")
}else{family_suffs = c("P")}	
						 
analnum = 1

#Get the information for the suspect intervals from the genomeCoverageBed files for both probands and parents
for (analysis_name in analyses){
	for (family_suff in family_suffs) {
		if (analtype == 'trio'){
			analysis = paste(analysis_name,family_suff,sep="_")
		}else{	analysis = analysis_name}	
		print(analysis)
		geno_cov_f <- paste(analysis,gen_cov_suffix,sep="_")
		#Check that the file exists and is non empty
		if (file.info(geno_cov_f)$size > 0 & file.exists(geno_cov_f)){
			geno_cov <- read.delim(geno_cov_f,stringsAsFactors=F,head=F,quote="",sep="\t")
			#Get a composite id with chr_start_end of the exon
			
			geno_cov$compid <- paste(geno_cov$V1,paste(geno_cov$V2,geno_cov$V3,sep="-"),sep=":")
			#Select columns to investigate
			geno_cov_filt <- geno_cov[,c("compid","V4")]
			colnames(geno_cov_filt) <- c("compid",analysis)
			#Define the complete file if it is empty
			if (analnum == 1){
				geno_cov_all <- geno_cov_filt
				
			}else{
				#Merge with the all file	
				geno_cov_all <- merge(geno_cov_all,geno_cov_filt,by="compid",all=T)
			}	
			analnum = analnum+1						
		}

	}
}

#Set all NA to 0
geno_cov_all[is.na(geno_cov_all)] <- 0
#geno_cov_all$compid<-NULL

geno_cov_all_n<-as.matrix(sapply(geno_cov_all[, 2:ncol(geno_cov_all)], as.numeric))
rownames(geno_cov_all_n)<-geno_cov_all$compid

het_table <- data.frame(
	 samplename=character(), 
	 compid=character(),
	  preads=character(), 
	  avgreads=character(),  
	 stringsAsFactors=FALSE) 
dups_table <- data.frame(
	 samplename=character(), 
	 compid=character(),
	  preads=character(), 
	  avgreads=character(),  
	 stringsAsFactors=FALSE) 
	 
for (interv_num in 1:nrow(geno_cov_all_n) ){
	#geno_cov_all_n_f <- as.data.frame(geno_cov_all_n[interv_num,geno_cov_all_n[interv_num,] >= outl_thr])
	upperbound <- mean(geno_cov_all_n[interv_num,])/2 + sd(geno_cov_all_n[interv_num,])/6
	lowerbound <- mean(geno_cov_all_n[interv_num,])/2 - sd(geno_cov_all_n[interv_num,])/6
	het_del <- geno_cov_all_n[interv_num,geno_cov_all_n[interv_num,] >= lowerbound & geno_cov_all_n[interv_num,]<= upperbound]
	doublemean <- mean(geno_cov_all_n[interv_num,])*2
	dups <- geno_cov_all_n[interv_num,geno_cov_all_n[interv_num,] >= doublemean]
	print((rownames(geno_cov_all_n))[interv_num])
	
	#For each interval check how many samples satisfies the condition of number of reads that are half
	#compared with the mean. If only one satisfies, then print:
	#interval, sample, number of reads, average number of reads
	#if (length(het_del) < 3 & length(het_del) > 0){
	if (length(het_del) > 0){
		#To get the sample name here I have to make two different procedures
		#because when there is a single result R will get a different data structure
		if (length(het_del) == 1){
			samplename = colnames(geno_cov_all_n)[which(geno_cov_all_n[interv_num,]==het_del[[1]])]
			preads <- het_del[[1]]
		}else{
			for (het_deln in 1:length(het_del)){
				samplename <- names(het_del)[het_deln]
				preads <- het_del[[het_deln]]			
			}
			#print (names(het_del))
			#print (interv_num)
		}
		#Get the compid
		compid <- rownames(geno_cov_all_n)[interv_num]
		#Get the average number of reads
		avgreads <- mean(geno_cov_all_n[interv_num,])
		#Fill the table with this result
		het_table[nrow(het_table) + 1,] = list(samplename,compid,preads,avgreads)	
		#selsample <- names(het_del)[1]
		#print(selsample)
	}	
	
	###Duplications
	if (length(dups) > 0){
		#To get the sample name here I have to make two different procedures
		#because when there is a single result R will get a different data structure
		if (length(dups) == 1){
			samplename = colnames(geno_cov_all_n)[which(geno_cov_all_n[interv_num,]==dups[[1]])]
			preads <- dups[[1]]
		}else{
			for (dupsn in 1:length(dups)){
				samplename <- names(dups)[dupsn]
				preads <- dups[[dupsn]]			
			}
			#print (names(dups))
			#print (interv_num)
		}
		#Get the compid
		compid <- rownames(geno_cov_all_n)[interv_num]
		#Get the average number of reads
		avgreads <- mean(geno_cov_all_n[interv_num,])
		#Fill the table with this result
		dups_table[nrow(dups_table) + 1,] = list(samplename,compid,preads,avgreads)	
		#selsample <- names(het_del)[1]
		#print(selsample)
	}	
}

max_cnv_num = 25
#Filter the HET TABLE removing all those samples having more than "max_cnv_num" cnvs
hd_freq_table <- as.data.frame(ftable(het_table$samplename))
hd_freq_table_f <- hd_freq_table[hd_freq_table$Freq<max_cnv_num,]
het_table_f <- merge(het_table,hd_freq_table_f,by.x="samplename",by.y="Var1",all.x=F,all.y=F)
	
#Filter the DUPS TABLE removing all those samples having more than "max_cnv_num" cnvs
dups_freq_table <- as.data.frame(ftable(dups_table$samplename))
dups_freq_table_f <- dups_freq_table[dups_freq_table$Freq<max_cnv_num,]
dups_table_f <- merge(dups_table,dups_freq_table_f,by.x="samplename",by.y="Var1",all.x=F,all.y=F)
