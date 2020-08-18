export function stringToUTCString(value: string): string {
  return new Date(value).toUTCString();
}

const roleMap = {
  100: "Admin",
  10: "Standard",
  0: "Guest",
};

export function roleToString(value: number): string {
  if (value === -1) {
    return "";
  }
  return roleMap[value];
}
