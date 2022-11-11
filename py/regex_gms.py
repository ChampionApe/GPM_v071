import re

def macros_from_str(string):
	a = re.findall(r"^.*?\$MACRO(.*?)\n",string,re.IGNORECASE | re.DOTALL | re.MULTILINE)
	b = ['$MACRO' + a[i] for i in range(len(a))]
	return {macro_name_from_str(b[i]): b[i] for i in range(len(b))}
	
def functions_from_str(string):
	a = re.findall(r"^.*?\$FUNCTION(.*?)\$ENDFUNCTION.*?$",string,re.IGNORECASE | re.DOTALL | re.MULTILINE)
	b = ['$FUNCTION' + a[i] + '$ENDFUNCTION' for i in range(len(a))]
	return {function_name_from_str(b[i]): b[i] for i in range(len(b))}
	
def blocks_from_str(string):
	a = re.findall(r"^.*?\$BLOCK(.*?)\$ENDBLOCK.*?$",string,re.IGNORECASE | re.DOTALL | re.MULTILINE)
	if a:
		return [block_name_from_str('$BLOCK'+a[i] +'$ENDBLOCK') for i in range(len(a))]
	else:
		return []
def groups_from_str(string):
	a = re.findall(r"^.*?\$GROUP(.*?)\;.*?$",string,re.IGNORECASE | re.DOTALL | re.MULTILINE)
	if a:
		return [group_name_from_str('$GROUP'+a[i] +';') for i in range(len(a))]
	else:
		return []

def macro_name_from_str(string):
	return re.search(r"^.*?\$MACRO(.*?)\(",string,re.IGNORECASE | re.DOTALL | re.MULTILINE).group(1).strip()
def function_name_from_str(string):
	return re.search(r"^.*?\$FUNCTION(.*?)\(.*?",string,re.IGNORECASE | re.DOTALL | re.MULTILINE).group(1).strip()
def block_name_from_str(string):
	return re.search(r"^.*?\$BLOCK\s(.*?)\s.*?",string,re.IGNORECASE | re.DOTALL | re.MULTILINE).group(1).strip()		
def group_name_from_str(string):
	return re.search(r"^.*?\$GROUP\s(.*?)\s.*?",string,re.IGNORECASE | re.DOTALL | re.MULTILINE).group(1).strip()		
