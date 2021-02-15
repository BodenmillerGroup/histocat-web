import {
  Alignment,
  Button,
  Menu,
  MenuItem,
  Navbar,
  NavbarDivider,
  NavbarGroup,
  NavbarHeading,
} from "@blueprintjs/core";
import { Popover2 } from "@blueprintjs/popover2";

import { useAuthStore } from "modules/auth";
import { appName } from "../env";
import { useHistory } from "react-router-dom";

export function MainNavBar() {
  const { userLogout } = useAuthStore();
  const history = useHistory();

  return (
    <Navbar fixedToTop={false}>
      <NavbarGroup align={Alignment.LEFT}>
        <NavbarHeading>{appName}</NavbarHeading>
        <NavbarDivider />
        <Button minimal={true} icon="home" text="Home" onClick={() => history.push("/")} />
        <Button minimal={true} icon="user" text="Users" onClick={() => history.push("/admin/users")} />
        <Button minimal={true} icon="predictive-analysis" text="Models" onClick={() => history.push("/admin/models")} />
      </NavbarGroup>
      <NavbarGroup align={Alignment.RIGHT}>
        <Popover2
          content={
            <Menu>
              <MenuItem text="Profile" href="/profile" />
              <MenuItem text="Logout" onClick={userLogout} />
            </Menu>
          }
          placement="bottom"
        >
          <Button minimal={true} icon="user" />
        </Popover2>
        <Button minimal={true} icon="notifications" />
        <Button minimal={true} icon="cog" />
      </NavbarGroup>
    </Navbar>
  );
}
