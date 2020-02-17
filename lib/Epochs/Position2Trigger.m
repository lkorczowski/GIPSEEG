function Trigger=Position2Trigger(StimPos,LengthTrigger)
%convert the input StimPos to Trigger
% StimPos contains only the position of the non-null element of Trigger
% Trigger is a channel containing a non-null element at the position of the
% stimulation

Trigger=zeros(LengthTrigger,1);
Trigger(StimPos)=1;
