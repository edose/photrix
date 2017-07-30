# Python Love/Hate

**Three things I love about python:***

 1. **Unit testing.** Python has it, R doesn't--simple as that. I thought my photometry workflow might stay simple enough to do without this. Ha. And now I can have real confidence in the 'simple' stuff on which my increasingly complex stuff is built.
 2. **Object-orientation / real classes**: Classes greatly help focus attention. They allow complex code structures without losing your way (read: mind). This is also sufficient reason, by itself, to port.
 3. **Execution speed**: Nice to have. R is blindingly fast for the vector stuff it's meant for, but pretty slow for general programming.
 4. **PyCharm pro IDE**: Luscious-looking. Sometimes clumsy with package management (site vs virtualenv), and nags you when *it* gets typing wrong. But wonderful unit-testing support, and did I mention that it Looks. So. Good.

OK, four--four things I love about python. 
And now for...

**Three things I dearly hate about python:**

 1. **No private variables, methods, functions**, etc: Oh, Jesus how could they have been so ****ing stupid? Stupid as the day is long. No more to say--just stupid.
Except that I will say more (of course). Why not just build your house in the middle of a freeway and let everyone else drive through it all day and night? The only justification for prohibiting python coders from having a tool as essential as private objects is: "Don't you trust your fellow coders?" To which I answer: "No, of course I don't--whyever would you?--and furthermore what business is it of yours if I don't?"
And did I mention: it's stupid.
 2. **Platform schizophrenia / coder abuse**: The open-source model is probably at fault here: everyone tries to inflict pet prejudices (read: errors) on the rest of the world. Example: pandas `df1.extend(df2)` returns a copy of the extended dataframe. Naturally. But python's own `list1.extend(list2)` returns...nothing. 
So: `df1.extend(df2).extend(df3)` does exactly what you'd expect. 
Whereas the natural syntax `list1.extend(list2).extend(list3)` stops your program cold. Well, excuse me for being reasonable. 
And don't even mention statsmodels trying to cram 'exog' and 'endog' variables down their users' throats, when at 2 standard naming systems have been accepted worldwide *for a century*. And statsmodels mixed-model predict() forgets to include random effects in their predictions at all, even as an option--and their docs *brag* about this failure.
Python is littered with this crap, everywhere, often costing hours to figure out and code around. They could learn a lot from R and even (horrors) from .Net.
 3. **'Self' is plastered everywhere**. This is so self-absorbed. Me, me, me. Self, self, self. So perfect for the age of Trump, so perfect for a 30-year-old language still living in its parent's basement.