const maybeTruncateString = (str: string, maxLength: number) => {
  let truncatedString: string | null = null;
  if (str.length > maxLength) {
    truncatedString = `${str.slice(0, maxLength / 2)}…${str.slice(-maxLength / 2)}`;
  }

  return truncatedString;
};

export default maybeTruncateString;
