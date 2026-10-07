function [OutputTemplate] = mc_GenPath_lite(Template,Subject,Task,Run,Session)

    OutputTemplate = strrep(Template,'[Subject]',Subject);
    OutputTemplate = strrep(OutputTemplate,'[Session]',Session);
    OutputTemplate = strrep(OutputTemplate,'[Run]',Run);
    OutputTemplate = strrep(OutputTemplate,'[Task]',Task);
    