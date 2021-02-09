%% Automatic documentation generation based on index.m file
% <index.html *INDEX*>
%
% INPUT
%
% OUTPUT
%
% DEPENDENCE: <index.html index>
%
% AUTHOR: J.-Y. Baudais, October 15, 2020

%% User guide
% *1.* Put this file on your working directory where |index.m| file is.
% |index.m| file contains only comments, it looks like a table of contents
% of the Matlab files.
%
% *2.* |index.m| file must contain the four lines
%
%     % *Documentation*
%     % 
%     % * <index.html index>
%     % * <generate_doc.html generate_doc>
%
% at then end.
%
% *3.* All lines in |index.html| with |<| are parsed to get the name of the
% file to be published. So, don't use |<| in line without filename to
% publish!
%
% *4.* The format of the parsed line must be (take care of all the *24*
% characters of this following line)
% 
%
%     % * <myfile.html myfile>
%
% where |myfile| has to be replaced by the name of your file. Write as many
% lines as files. Comments can be added in |index.m| using sectionning to
% organized the index file. Don't use space, |<| nor |.html| in |myfile|
% name. The skeleton of the header of |myfile.m| is
%
%     %% One line description of the code
%     % <index.html *INDEX*>
%     % 
%     % Long description of the code with references and equations
%     % 
%     % INPUT
%     % 
%     %     X: description
%     %     Y: description 
%     %     ...
%     % 
%     % OUTPUT
%     % 
%     %     Z: description
%     %     T: description
%     %     ...
%     % 
%     % DEPENDENCE: <other_file.html other_file.m> ; <another_one.html
%     % another_one.m>...
%     % 
%     % AUTHOR: name and date
%
% INPUT, OUTPUT or DEPENDENCE can be empty (see the header of this
% |generate_doc.m| file).
%
% *5.* See "Publishing Markup" documentation to write comments.
%
% *6.* Check your all your |.m| files in the current folder, click on
% "Sow Current Folder Actions" in the current folder browser and select
%
%    1. Reports > Code Analyzer Report
%    2. Reports > Dependency Report
%
% *7.* Run this |generate_doc.m| file.
%
% *Do not change the code after this line without knowing what you are
% doing*

%% The code

clear;
close all;
clc;

fileName='index.m';
fId = fopen(fileName);
y = 0;
tLine = fgetl(fId);
while ischar(tLine)
   num = length(strfind(tLine, '>'));
   if num > 0
      y = y + num;
      tmpName=textscan(tLine, '%% %s %s %s');
      %fPublish=regexprep(char(tmpName{3}),'>','');
      fPublish=regexprep(regexprep(char(tmpName{2}),'<',''),'.html','');
      if fopen(strcat(fPublish,'.m'))==-1
          fprintf('%s.m: not exist <===== check file name in index.m\n',fPublish);
      else
          publish(fPublish,struct('outputDir','./html','evalCode',false));
          fprintf(1,'%s: published\n',fPublish);
      end
   end
   tLine = fgetl(fId);
end
fclose(fId);
web('html/index.html')
