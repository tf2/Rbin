`RbinConvert_aCGH` <- function(inputname, binname) {
	.C("makeaCGHbin", "inputname"=as.character(inputname), "biname"=as.character(binname))
}

`RbinConvert_exome` <- function(inputname, binname) {
	.C("makeexomebin", "inputname"=as.character(inputname), "biname"=as.character(binname))
}

`RbinRead_aCGH` <- function(chr, start, stop, indexfile = NULL, filenames=NULL, verbose=TRUE) {
	
	offsets <- .C( "getOffsets",
			"filename" = as.character(indexfile),
			"chr" = as.integer(chr), 
	  		"start" =as.integer(start),
	   		"stop" = as.integer(stop),
	    	"rows" = as.integer(0),
	    	"offset" = as.double(0)
			)
	if(verbose) {
		print("Calculated Offsets and Length!!!")
	}
	positions <- .C( "getPositions",
			"filename" = as.character(indexfile),
			"chr" = as.integer(chr), 
	  		"start" =as.integer(start),
	   		"stop" = as.integer(stop),
	    	"chr_values" = as.integer(1:offsets$rows),
	    	"start_values" = as.integer(1:offsets$rows),
	    	"stop_values" = as.integer(1:offsets$rows)
			)
	if(verbose) {
		print("Retrived Genomic Positions!!!")
	}
	matr = matrix(nrow=offsets$rows, ncol=length(filenames))
	if(verbose) {
		pb <- txtProgressBar(min = 0, max = length(filenames), style = 3);
	}
	for(x in 1:length(filenames)) {
		rr= .C( "getValues",
			"filename" = as.character(filenames[x]),
	 		"start_position" = as.double(offsets$offset), 
	  		"number_row_to_read" =as.integer(offsets$rows),
	   		"values" = as.double(1:offsets$rows)
			)
			matr[,x] = rr$values; 
			if(verbose) {
				setTxtProgressBar(pb, x); rr=NULL;
			}
	   } 
	   data = cbind(positions$chr_values, positions$start_values, positions$stop_values, matr)
	   colnames(data)=c("chr", "start", "stop", filenames)
		if(verbose) {   
			close(pb);
		}
	if(verbose) {
		print(paste("Retirved ", dim(data)[1]*dim(data)[2], " datapoints with ", length(filenames), " sequential file reads!!!",sep=""))
	}
	
	positions=NULL; offsets=NULL; rr=NULL;
	invisible(data)
}

`RbinRead_exome` <- function(chr, start, stop, typ="adm3score", indexfile = NULL, filenames=NULL) {
	
	offsets <- .C( "geteOffsets",
				  "filename" = as.character(indexfile),
				  "chr" = as.integer(chr), 
				  "start" =as.integer(start),
				  "stop" = as.integer(stop),
				  "rows" = as.integer(0),
				  "offset" = as.double(0)
				  )
	print("Calculated Offsets and Length!!!")
	
	positions <- .C( "getePositions",
					"filename" = as.character(indexfile),
					"chr" = as.integer(chr), 
					"start" =as.integer(start),
					"stop" = as.integer(stop),
					"chr_values" = as.integer(1:offsets$rows),
					"start_values" = as.integer(1:offsets$rows),
					"stop_values" = as.integer(1:offsets$rows)
					)
	print("Retrived Genomic Positions!!!")
	
	matr = matrix(nrow=offsets$rows, ncol=length(filenames))
	pb <- txtProgressBar(min = 0, max = length(filenames), style = 3);
	
	for(x in 1:length(filenames)) {
		rr = NULL
		if(typ=="adm3score") {
			rr= .C( "geteValues2",
			   "filename" = as.character(filenames[x]),
			   "start_position" = as.double(offsets$offset), 
			   "number_row_to_read" =as.integer(offsets$rows),
			   "values" = as.double(1:offsets$rows)
			   )
			
		} else {
			rr= .C( "geteValues1",
				   "filename" = as.character(filenames[x]),
				   "start_position" = as.double(offsets$offset), 
				   "number_row_to_read" =as.integer(offsets$rows),
				   "values" = as.double(1:offsets$rows)
				   )
			
		}
		matr[,x] = rr$values; setTxtProgressBar(pb, x); rr=NULL;
	} 
	data = cbind(positions$chr_values, positions$start_values, positions$stop_values, matr)
	colnames(data)=c("chr", "start", "stop", filenames)
	close(pb);
	print(paste("Retirved ", dim(data)[1]*dim(data)[2], " datapoints with ", length(filenames), " sequential file reads!!!",sep=""))
	positions=NULL; offsets=NULL; rr=NULL;
	invisible(data)
}

`Rbin_example` <- function(chr=1, start=1, stop=20000000, n=1000) {
	
	# Make the example binary file format
	path = system.file("extdata", package="Rbin");
	binname = paste(path, "/", "bin.dat", sep=""); inname = paste(path, "/", "Example_Input.txt", sep="");
	RbinConvert_aCGH(inname,binname);
	
	# Make the example index file and ensure (just for example) that it is sorted correctly.
	index = read.table(inname); index=index[order(index[,1], index[,2], index[,3]),];
	indexfile = paste(path, "/", "Example_Index.txt", sep="");
	write.table(index[,1:3], indexfile, sep="\t", row.names=F, col.names=F, quote=F);
	
	# Replicate the binary filename by n and execute n sequetial file reads returning a matrix for the given position.
	filenames = rep(binname, n);
	r = RbinRead_aCGH(chr, start, stop, indexfile, filenames);
	
	invisible(r)
}
