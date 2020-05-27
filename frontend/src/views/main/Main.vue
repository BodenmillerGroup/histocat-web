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
            <v-list-item to="/main/experiments">
              <v-list-item-icon>
                <v-icon>mdi-view-dashboard-outline</v-icon>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Dashboard</v-list-item-title>
              </v-list-item-content>
            </v-list-item>
          </v-list>
          <v-divider />
          <v-list nav dense v-if="activeExperimentId">
            <v-list-item @click="showData(false)">
              <v-list-item-icon>
                <v-icon>mdi-magnify-scan</v-icon>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Image</v-list-item-title>
              </v-list-item-content>
            </v-list-item>
            <v-list-item @click="showData(true)">
              <v-list-item-icon>
                <v-icon>mdi-scatter-plot-outline</v-icon>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Data</v-list-item-title>
              </v-list-item-content>
            </v-list-item>
          </v-list>
          <v-divider v-if="isAdmin && activeExperimentId" />
          <v-list nav dense v-if="isAdmin">
            <v-list-item to="/main/admin/users/all">
              <v-list-item-icon>
                <v-icon>mdi-account-multiple-outline</v-icon>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Users</v-list-item-title>
              </v-list-item-content>
            </v-list-item>
            <v-list-item to="/main/admin/experiments/all">
              <v-list-item-icon>
                <v-icon>mdi-folder-multiple-outline</v-icon>
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Experiments</v-list-item-title>
              </v-list-item-content>
            </v-list-item>
          </v-list>
          <v-list nav dense flat>
            <v-list-item @click="switchMiniDrawer">
              <v-list-item-icon>
                <v-icon v-html="miniDrawer ? 'mdi-chevron-right' : 'mdi-chevron-left'" />
              </v-list-item-icon>
              <v-list-item-content>
                <v-list-item-title>Collapse</v-list-item-title>
              </v-list-item-content>
            </v-list-item>
          </v-list>
        </v-col>
      </v-row>
    </v-navigation-drawer>
    <v-app-bar app dense dark color="primary" :clipped-left="$vuetify.breakpoint.lgAndUp" extension-height="0">
      <v-app-bar-nav-icon @click.stop="switchShowDrawer" />
      <v-toolbar-title @click="$router.push('/main/experiments')" class="toolbar-title">{{ appName }}</v-toolbar-title>
      <v-spacer />
      <v-btn-toggle v-model="views" multiple background-color="primary" group>
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
    <v-content>
      <router-view />
    </v-content>
  </div>
</template>

<script lang="ts">
import ToolbarProgressBar from "@/components/ToolbarProgressBar.vue";
import { appName } from "@/env";
import { mainModule } from "@/modules/main";
import { BroadcastManager } from "@/utils/BroadcastManager";
import { WebSocketManager } from "@/utils/WebSocketManager";
import { Component, Vue, Watch } from "vue-property-decorator";
import { experimentModule } from "@/modules/experiment";

const routeGuardMain = async (to, from, next) => {
  if (to.path === "/main") {
    next("/main/experiments");
  } else {
    next();
  }
};

@Component({
  components: { ToolbarProgressBar },
})
export default class Main extends Vue {
  readonly experimentContext = experimentModule.context(this.$store);
  readonly mainContext = mainModule.context(this.$store);

  appName = appName;
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

  get activeExperimentId() {
    return this.experimentContext.getters.activeExperimentId;
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

  showData(value: boolean) {
    this.mainContext.mutations.setShowData(value);
  }

  switchShowDrawer() {
    this.mainContext.mutations.setDashboardShowDrawer(!this.mainContext.getters.dashboardShowDrawer);
  }

  switchMiniDrawer() {
    this.mainContext.mutations.setDashboardMiniDrawer(!this.mainContext.getters.dashboardMiniDrawer);
  }

  get isAdmin() {
    return this.mainContext.getters.hasAdminAccess;
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
.subheader {
  font-size: 10px;
  font-weight: bold;
}

.toolbar-title {
  cursor: pointer;
}
</style>
