%LET PATH1=/folders/myfolders/relvm/output4133;
LIBNAME HC "&PATH1";

FILENAME REFFILE '/folders/myfolders/relvm/input/dat4133.csv';

PROC IMPORT DATAFILE=REFFILE
	DBMS=CSV replace
	OUT=HC.dat;
	GETNAMES=YES;
RUN;

proc nlmixed data=HC.dat tech=dbldog qpoints=30 noad;
	parms mu1-mu7=0, gamma1-gamma7=1, err1-err7=1;
	array x[7] x1-x7;
	array mu[7] mu1-mu7;
	array gamma[7] gamma1-gamma7;
	array log_den[7] log_den1-log_den7;
	array w[7] w1-w7;
	array err[7] err1-err7;
	logmu=0;

	do i=1 to 7;

		if x{i}=. then
			log_den{i}=0;
		else
			log_den{i}=w{i}*log(pdf('normal', x{i}, mu{i}+gamma{i}*alpha, err{i}));
		logmu=log_den{i} + logmu;
	end;
	model pid ~ general(logmu);
	random alpha ~ normal(0, 1) subject=pid;
	ods output ParameterEstimates=HC.Out_Est;
	predict alpha out=HC.Out_Pred;
run;