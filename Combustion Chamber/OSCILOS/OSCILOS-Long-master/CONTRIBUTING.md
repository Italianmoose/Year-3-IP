OSCILOS is being developed for research. Please [join us](http://www.oscilos.com/)!. 

## Git and Pull requests
* Contributions are welcome. For outside users, they can be submitted with GitHub pull request, and will be reviewed and accepted by the team. The project uses the Fork & Pull model; for more details see [here](https://help.github.com/articles/using-pull-requests).
* Internal developers have write access to the branches, and can merge changes themselves. They do not need to make Pull requests.
* The project uses the [GitFlow branching model](http://nvie.com/posts/a-successful-git-branching-model/). As such please make a new branch for every feature you are working on.
* Internal developers should merge Bugfixes and Hotfixes to `master-private` and `develop`.
* After a release is made to `master-private`, or a Bugfix/Hotfix is implemented, changes in `master-private` should be merged in the public repo. `master` branch.
* The latest changes available to users outside of MorgansLab are in the public repo. `master` branch. To have access to newer bug fixes or become a collaborator, contact the repo. admin. 
* Make sure your commit messages are clear. The first line should be a summary of the commit, following lines should contain details about the mentioned changes. 
* Please avoid using re-basing and fast forward merges. 

## Coding style and version numbering
* Use a lot of comments; more is better than not enough. 
* Indent your code using the Matlab automated indent feature. 
* When incrementing the version number, please stick to the [Semantic Versioning Standards](http://semver.org/)

## Documentation
* The documentation can be found in the [docs](docs) folder. 

## Logging issues and fixing issues
* Please log any issues you find [here](https://github.com/MorgansLab/OSCILOS-1.2/issues), and don't forget to say what version of Matlab you are using, and what OS. 
* Label your issue as a `bug`, `enhancement` etc.
* Developers working on an issue should assign themselves the issue before starting any work. 
* Corrections for `bug` issues should be done using GitFlow Hotfix or Bugfix, or in a release branch. 
* When committing to correct an issue, use the wording `Fixes #1` or `Resolves #1` or `Closes #1`. Doing this will automatically [close issue Number 1, and link to the relevant commit](https://help.github.com/articles/closing-issues-via-commit-messages/)
* Use the template below when logging issues:

````
*Matlab Version : r2014a*
*OS : Windows 8.1*
*Commit: 86f0a46d3e01461926622151eb9fff689c677321*

## Observed behaviour
When selecting a temperature ratio across the flame of 1 (not a plausible configuration, but hey) the plotting shows error

### Error message

\`\`\`\`
Error using set
Bad property value found.
Object Name: axes
Property Name: 'YLim'
Value must be a 2 element vector.

Error in GUI_INI_TP>Fcn_GUI_INI_TP_PLOT (line 1250)
        set(hAxes1,'ylim',[ymin ymax])

Error in GUI_INI_TP>pb_Plot_Callback (line 480)
Fcn_GUI_INI_TP_PLOT(hObject)

Error in gui_mainfcn (line 95)
        feval(varargin{:});

Error in GUI_INI_TP (line 42)
    gui_mainfcn(gui_State, varargin{:});

Error in
@(hObject,eventdata)GUI_INI_TP('pb_Plot_Callback',hObject,eventdata,guidata(hObject))

 
Error using waitfor
Error while evaluating uicontrol Callback
\`\`\`\`

## Expected behaviour
Ylim should be adjusted when the plot is a horizontal line
````