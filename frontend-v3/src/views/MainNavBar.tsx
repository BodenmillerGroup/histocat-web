import { Alignment, Button, Navbar, NavbarDivider, NavbarGroup, NavbarHeading } from "@blueprintjs/core";

export function MainNavBar() {
  return (
    <Navbar fixedToTop={true}>
      <NavbarGroup align={Alignment.LEFT}>
        <NavbarHeading>histoCAT</NavbarHeading>
        <NavbarDivider />
        <Button minimal={true} icon="home" text="Home" />
        <Button minimal={true} icon="document" text="Files" />
      </NavbarGroup>
      <NavbarGroup align={Alignment.RIGHT}>
        <Button minimal={true} icon="user" />
        <Button minimal={true} icon="notifications" />
        <Button minimal={true} icon="cog" />
      </NavbarGroup>
    </Navbar>
  );
}
