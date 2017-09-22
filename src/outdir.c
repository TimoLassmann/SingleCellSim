
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "tldevel.h"
#include "outdir.h"

#include <sys/types.h>
#include <sys/stat.h>

#include <dirent.h>
#include <unistd.h>

char* concat_out_dirs(char* root, char* sub)
{
	char* tmp = NULL;
	int len = (int)(strlen(root) + strlen(sub) + 2 );
	MMALLOC(tmp, sizeof(char) *len);
	snprintf(tmp,len,"%s/%s",root, sub);
	
	
	return tmp;
ERROR:
	return NULL;
}

int create_dir(char* name, int force)
{
	struct stat st = {0};
	
	if(stat(name, &st) == -1) {
		RUN(mkdir(name, 0700));
	}else{
		if(!force){
			ERROR_MSG("Directory %s already exists.",name);
		}
	}
	return OK;
ERROR:
	return FAIL;
}

int set_log_file(char* root, char* name)
{
	char buffer[BUFFER_LEN];
	char out_root_dir[BUFFER_LEN];
	snprintf(out_root_dir, BUFFER_LEN, "%s/" ,root);
	snprintf(buffer, BUFFER_LEN, "%s%s/%s.log",out_root_dir,OUTDIR_LOG,name);
	
	FILE* test_ptr = NULL;
	
	RUNP(test_ptr = fopen(buffer, "w"));
	fclose(test_ptr);
	
	tlog.set_logfile(buffer);
	
	
	
	return OK;
ERROR:
	WARNING_MSG("Log file: %s cannot be opened",buffer);
	return FAIL;
}

int create_output_directories(char* root)
{
	//struct stat st = {0};
	char buffer[BUFFER_LEN];
	char out_root_dir[BUFFER_LEN];
	
	//create main output directory.
	
	snprintf(out_root_dir, BUFFER_LEN, "%s/" ,root);
	
	RUN(create_dir(out_root_dir,1));
	
	/*if(stat(out_root_dir, &st) == -1) {
		RUN(mkdir(out_root_dir, 0700));
	}else{
		ERROR_MSG("Directory %s already exists.",out_root_dir);
	}*/
	
	snprintf(buffer, BUFFER_LEN, "%s%s/",out_root_dir,OUTDIR_PARAM);
	RUN(create_dir(buffer,1));
        	
	snprintf(buffer, BUFFER_LEN, "%s%s/",out_root_dir,OUTDIR_CHECKPOINTS);
	RUN(create_dir(buffer,1));
	
	snprintf(buffer, BUFFER_LEN, "%s%s/",out_root_dir,OUTDIR_LOG);
	RUN(create_dir(buffer,1));

	
		
	
	// create sub_directories;
	
	/*for(i =0; i < 100;i++){
		snprintf(buffer, BUFFER_LEN, "%s%s/%02d/",out_root_dir,OUTDIR_RAWSEQ, i  );
		//fprintf(stdout,"I want to create directory: %s\n", buffer);
		if (stat(buffer, &st) == -1) {
	 
	 RUN(mkdir(buffer, 0700));
		}
		for(j =0; j < 100;j++){
	 snprintf(buffer, BUFFER_LEN, "%s%s/%02d/%02d/",out_root_dir,OUTDIR_RAWSEQ, i,j  );
		//	fprintf(stdout,"I want to create directory: %s\n", buffer);
	 if (stat(buffer, &st) == -1) {
	 RUN(mkdir(buffer, 0700));
	 }
		}
	 }*/
	
	return OK;
ERROR:
	return FAIL;
}

int check_if_output_directories_exists(char* root)
{
	struct stat st = {0};
	char buffer[BUFFER_LEN];
	char out_root_dir[BUFFER_LEN];
	
	snprintf(out_root_dir, BUFFER_LEN, "%s/" ,root);
	//test if root is there
	if (stat(out_root_dir, &st) != -1) {	
		snprintf(buffer, BUFFER_LEN, "%s%s/",out_root_dir,OUTDIR_CHECKPOINTS);
		if (stat(buffer, &st) == -1) {
			ERROR_MSG("%d not found.",buffer );
		}
	}else{
		return FAIL;
	}
	return OK;
ERROR:
	return FAIL;
}


// function to retrive file list ending in specific suffix.
// need a pointer to a char** ... and a num_item variabel...

char**  get_files_names(char* directory, char* suffix, char** list, int* num_files,int depth)
{
	DIR *dp;
	struct dirent *entry;
	struct stat statbuf;

	char cwd[BUFFER_LEN];
	
//ASSERT(list != NULL, " List is NULL");
	ASSERT(directory !=NULL," directory is NULL");
	ASSERT(suffix != NULL, "suffix is NULL");
	ASSERT(strlen(suffix) >1, "suffic is too short: %s", suffix);
	
	// get current directory...
	RUNP(getcwd(cwd, sizeof(cwd)));
	//open directory in question to count number of files...
	RUNP(dp = opendir(directory));
	
	if(chdir(directory)){
		ERROR_MSG("chdir failed.");
	}
	
	while((entry = readdir(dp)) != NULL) {
		lstat(entry->d_name,&statbuf);
		if(S_ISDIR(statbuf.st_mode)) {
			/* Found a directory, but ignore . and .. */
			if(strcmp(".",entry->d_name) == 0 ||
			   strcmp("..",entry->d_name) == 0)
				continue;
			//printf("%*s%s/\n",depth,"",entry->d_name);
			/* Recurse at a new indent level */
			/*
			strcpy(org_dir,directory);
			
			MREALLOC(directory, sizeof(char) * (strlen( directory) + strlen(entry->d_name) +2));
			strcat(directory, "/");
			strcat(directory, entry->d_name);
			fprintf(stdout,"TEST: %s org\n",org_dir);
			fprintf(stdout,"TEST: %s org\n",directory);
			
			list = get_files_names(directory,suffix,list,num_files,  depth+4);
			//turn directoty back to original
			directory[0] = 0;
			
			strcpy(directory,org_dir);
			
			fprintf(stdout,"TEST: %s back to org\n",directory);*/
			
		}else{
			if(!strcmp(suffix, entry->d_name+ (strlen(entry->d_name) - strlen(suffix)))){
				
				//fprintf(stdout,"REALLOCING: %d size:%d\n", *num_files, (int)( sizeof(char*) * ( *num_files)));
				MREALLOC(list, sizeof(char*) * ( *num_files + 1));
				//fprintf(stdout,"%p List\n", list);
				list[ *num_files] = NULL;
				MMALLOC(list[ *num_files], sizeof(char) *(strlen(entry->d_name)+1));
				strcpy(list[ *num_files], entry->d_name);
				//fprintf(stdout,"%p List\n", list);
				
				//for(j = 0; j <  *num_files+1;j++){
				//	fprintf(stdout,"%d\t%s\n",j,list[j]);
				//}
				*num_files = *num_files + 1;
			}
		}
	}
	if(closedir(dp)){
		ERROR_MSG("closedir failed.");
	}
	if(chdir(cwd)){
		ERROR_MSG("chdir failed.");
	}
	return list;
ERROR:
	return NULL;
}

#ifdef ITEST

int main (int argc,char * argv[])
{
	
	char* dir = NULL;
	char** list = NULL;
	int num_items = 0;
	int i;
	
	if(argc != 3){
		fprintf(stdout,"Usage: ./outdir_ITEST <dir> <suffix>\n");
		return EXIT_SUCCESS;
	}
	
	MMALLOC(dir,sizeof(char) * strlen(argv[1])+1);
	strcpy(dir, argv[1]);
	
	list = get_files_names(dir, argv[2], list,&num_items, 0);
	if(num_items){
	
		fprintf(stdout,"Found %d %s files:\n\n",num_items,argv[2]);
	
		for(i = 0; i < num_items;i++){
			fprintf(stdout,"   %s\n",list[i]);
		}
	}else{
		ERROR_MSG("No files ending in %s found.\n",argv[2]);
	}
	
	
	for(i = 0; i < num_items;i++){
		MFREE(list[i]);
	}
	MFREE(dir);
	MFREE(list);
	
	return EXIT_SUCCESS;
ERROR:
	if(dir){
		MFREE(dir);
	}
	if(list){
		
		MFREE(list);
	}
	return EXIT_FAILURE;
}

#endif










