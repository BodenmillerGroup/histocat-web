<template>
  <div>
    <v-navigation-drawer
      :mini-variant="miniDrawer"
      mini-variant-width="60"
      :clipped="$vuetify.breakpoint.lgAndUp"
      v-model="showDrawer"
      fixed
      app
      width="180"
    >
      <v-row no-gutters>
        <v-col>
          <v-list nav dense>
            <v-list-item
              v-if="activeGroupId && !activeProjectId"
              :to="{ name: 'group-projects', params: { groupId: activeGroupId } }"
            >
              <v-list-item-icon>
                <v-tooltip right>
                  <template v-slot:activator="{ on }">
                    <v-icon v-on="on">mdi-view-dashboard-outline</v-icon>
                  </template>
                  <span>Projects</span>
                </v-tooltip>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Projects</v-list-item-title>
              </v-list-item-content>
            </v-list-item>

            <v-list-item
              v-if="activeGroupId && !activeProjectId"
              :to="{ name: 'group-members', params: { groupId: activeGroupId } }"
            >
              <v-list-item-icon>
                <v-tooltip right>
                  <template v-slot:activator="{ on }">
                    <v-icon v-on="on">mdi-account-multiple-outline</v-icon>
                  </template>
                  <span>Members</span>
                </v-tooltip>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Members</v-list-item-title>
              </v-list-item-content>
            </v-list-item>

            <v-list-item v-if="activeGroupId && activeProjectId" @click="showData(false)">
              <v-list-item-icon>
                <v-tooltip right>
                  <template v-slot:activator="{ on }">
                    <v-icon v-on="on">mdi-magnify-scan</v-icon>
                  </template>
                  <span>Image</span>
                </v-tooltip>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Image</v-list-item-title>
              </v-list-item-content>
            </v-list-item>

            <v-list-item v-if="activeGroupId && activeProjectId" @click="showData(true)">
              <v-list-item-icon>
                <v-tooltip right>
                  <template v-slot:activator="{ on }">
                    <v-icon v-on="on">mdi-scatter-plot-outline</v-icon>
                  </template>
                  <span>Data</span>
                </v-tooltip>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Data</v-list-item-title>
              </v-list-item-content>
            </v-list-item>

            <v-divider v-if="activeGroupId && activeProjectId" />

            <v-list-item
              v-if="activeGroupId && activeProjectId"
              @click="$router.push({ name: 'group-projects', params: { groupId: activeGroupId } }, () => {})"
            >
              <v-list-item-icon>
                <v-tooltip right>
                  <template v-slot:activator="{ on }">
                    <v-icon v-on="on">mdi-view-dashboard-outline</v-icon>
                  </template>
                  <span>Projects</span>
                </v-tooltip>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Projects</v-list-item-title>
              </v-list-item-content>
            </v-list-item>

            <v-list-item v-if="!activeGroupId && !activeProjectId && isAdmin" :to="{ name: 'admin-users' }">
              <v-list-item-icon>
                <v-tooltip right>
                  <template v-slot:activator="{ on }">
                    <v-icon v-on="on">mdi-account-outline</v-icon>
                  </template>
                  <span>Users</span>
                </v-tooltip>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Users</v-list-item-title>
              </v-list-item-content>
            </v-list-item>
          </v-list>
        </v-col>
      </v-row>
    </v-navigation-drawer>
    <v-app-bar app dense dark color="primary" :clipped-left="$vuetify.breakpoint.lgAndUp" extension-height="0">
      <v-app-bar-nav-icon @click.stop="switchShowDrawer" />
      <v-toolbar-title @click.stop="$router.push({ name: 'groups' }, () => {})" class="toolbar-title">{{
        appName
      }}</v-toolbar-title>
      <v-spacer />
      <v-btn-toggle v-if="activeGroupId && activeProjectId" v-model="views" multiple background-color="primary" group>
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" value="workspace" color="primary">
              <v-icon>mdi-file-tree</v-icon>
            </v-btn>
          </template>
          <span v-if="!showWorkspace">Show workspace</span>
          <span v-else>Hide workspace</span>
        </v-tooltip>
        <v-tooltip bottom>
          <template v-slot:activator="{ on }">
            <v-btn v-on="on" value="options" color="primary">
              <v-icon>mdi-tune</v-icon>
            </v-btn>
          </template>
          <span v-if="!showOptions">Show options</span>
          <span v-else>Hide options</span>
        </v-tooltip>
      </v-btn-toggle>
      <v-menu bottom left offset-y>
        <template v-slot:activator="{ on }">
          <v-btn icon v-on="on">
            <v-icon>mdi-dots-vertical</v-icon>
          </v-btn>
        </template>
        <v-list>
          <v-list-item to="/main/profile">
            <v-list-item-title>Profile</v-list-item-title>
            <v-list-item-action>
              <v-icon>mdi-account</v-icon>
            </v-list-item-action>
          </v-list-item>
          <v-list-item @click="logout">
            <v-list-item-title>Logout</v-list-item-title>
            <v-list-item-action>
              <v-icon>mdi-logout-variant</v-icon>
            </v-list-item-action>
          </v-list-item>
        </v-list>
      </v-menu>
      <ToolbarProgressBar
        :processing="processing"
        :progress="processingProgress"
        :indeterminate="false"
        color="light-blue lighten-2"
        slot="extension"
      />
    </v-app-bar>
    <v-main>
      <router-view />
    </v-main>
  </div>
</template>

<script lang="ts">
import ToolbarProgressBar from "@/components/ToolbarProgressBar.vue";
import { appName } from "@/env";
import { mainModule } from "@/modules/main";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { WebSocketManager } from "@/utils/WebSocketManager";
import { Component, Vue, Watch } from "vue-property-decorator";
import { groupModule } from "@/modules/group";
import { projectsModule } from "@/modules/projects";

const routeGuardMain = async (to, from, next) => {
  if (to.path === "/main") {
    next("/main/groups");
  } else {
    next();
  }
};

@Component({
  components: { ToolbarProgressBar },
})
export default class Main extends Vue {
  readonly mainContext = mainModule.context(this.$store);
  readonly groupContext = groupModule.context(this.$store);
  readonly projectsContext = projectsModule.context(this.$store);

  readonly appName = appName;
  views: string[] = ["workspace", "options"];

  beforeRouteEnter(to, from, next) {
    routeGuardMain(to, from, next);
  }

  beforeRouteUpdate(to, from, next) {
    routeGuardMain(to, from, next);
  }

  @Watch("views")
  viewsChanged(views: string[]) {
    this.mainContext.mutations.setLayout({
      showWorkspace: views.includes("workspace"),
      showOptions: views.includes("options"),
    });
  }

  get showWorkspace() {
    return this.mainContext.getters.showWorkspace;
  }

  get showOptions() {
    return this.mainContext.getters.showOptions;
  }

  get miniDrawer() {
    return this.mainContext.getters.dashboardMiniDrawer;
  }

  get showDrawer() {
    return this.mainContext.getters.dashboardShowDrawer;
  }

  set showDrawer(value: boolean) {
    this.mainContext.mutations.setDashboardShowDrawer(value);
  }

  switchShowDrawer() {
    this.mainContext.mutations.setDashboardShowDrawer(!this.mainContext.getters.dashboardShowDrawer);
  }

  get isAdmin() {
    return this.mainContext.getters.isAdmin;
  }

  get isGroupAdmin() {
    return this.groupContext.getters.isGroupAdmin;
  }

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  get activeProjectId() {
    return this.projectsContext.getters.activeProjectId;
  }

  async logout() {
    await this.mainContext.actions.userLogOut();
  }

  get processing() {
    return this.mainContext.getters.processing;
  }

  get processingProgress() {
    return this.mainContext.getters.processingProgress;
  }

  showData(value: boolean) {
    this.mainContext.mutations.setShowData(value);
  }

  mounted() {
    WebSocketManager.init(this.$store);
    BroadcastManager.init(this.$store);
  }

  beforeDestroy() {
    WebSocketManager.close();
    BroadcastManager.clear();
  }
}
</script>

<style scoped>
.toolbar-title {
  cursor: pointer;
}
</style>

<style>
.link {
  text-decoration: none;
}
</style>
