<template>
  <v-col>
    <v-toolbar dense class="toolbar">
      <v-toolbar-title>Projects</v-toolbar-title>
      <v-spacer />
      <v-toolbar-items>
        <v-btn
          v-if="isGroupAdmin"
          text
          color="primary"
          :to="{ name: 'group-projects-create', params: { groupId: activeGroupId } }"
        >
          Create Project
        </v-btn>
      </v-toolbar-items>
    </v-toolbar>
    <v-expansion-panels>
      <v-expansion-panel>
        <v-expansion-panel-header>Filter</v-expansion-panel-header>
        <v-expansion-panel-content>
          <v-select
            v-model="tagsFilter"
            :items="tags"
            item-text="text"
            item-value="value"
            chips
            clearable
            label="Tags"
            multiple
            prepend-icon="mdi-filter-outline"
            solo
            dense
          >
            <template v-slot:selection="{ attrs, item, select, selected }">
              <v-chip v-bind="attrs" :input-value="selected" close @click="select" @click:close="removeTagFilter(item)">
                {{ item }}
              </v-chip>
            </template>
          </v-select>
        </v-expansion-panel-content>
      </v-expansion-panel>
    </v-expansion-panels>
    <v-card>
      <v-card-title>
        <v-spacer />
        <v-text-field v-model="search" append-icon="mdi-magnify" label="Search" single-line hide-details clearable />
      </v-card-title>
      <v-data-table
        :headers="headers"
        :items="items"
        :loading="!items"
        :search="search"
        :custom-filter="filter"
        :items-per-page="15"
        :footer-props="{
          itemsPerPageOptions: [10, 15, 20, -1],
          showFirstLastPage: true,
          showCurrentPage: true,
        }"
        multi-sort
      >
        <template v-slot:item.name="{ item }">
          <router-link class="link" :to="{ name: 'group-project', params: { projectId: item.id } }">
            {{ item.name }}
          </router-link>
        </template>
        <template v-slot:item.tags="{ item }">
          <v-chip v-for="tag in item.tags" :key="tag" x-small pill disabled class="mr-1">
            {{ tag }}
          </v-chip>
        </template>
        <template v-slot:item.created_at="{ item }">
          {{ item.created_at | stringToUTCString }}
        </template>
        <template v-slot:item.action="{ item }">
          <v-menu bottom left>
            <template v-slot:activator="{ on }">
              <v-btn icon v-on="on">
                <v-icon>mdi-dots-vertical</v-icon>
              </v-btn>
            </template>
            <v-list dense>
              <v-list-item :to="{ name: 'group-projects-edit', params: { projectId: item.id } }">
                <v-list-item-icon>
                  <v-icon color="grey">mdi-pencil-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Edit</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
              <v-list-item v-if="isGroupAdmin" @click="deleteProject(item.id)">
                <v-list-item-icon>
                  <v-icon color="red accent-1">mdi-delete-outline</v-icon>
                </v-list-item-icon>
                <v-list-item-content>
                  <v-list-item-title>Delete</v-list-item-title>
                </v-list-item-content>
              </v-list-item>
            </v-list>
          </v-menu>
        </template>
      </v-data-table>
    </v-card>
  </v-col>
</template>

<script lang="ts">
import { projectsModule } from "@/modules/projects";
import { mainModule } from "@/modules/main";
import { Component, Vue } from "vue-property-decorator";
import { groupModule } from "@/modules/group";

@Component
export default class ProjectsView extends Vue {
  mainContext = mainModule.context(this.$store);
  groupContext = groupModule.context(this.$store);
  projectsContext = projectsModule.context(this.$store);

  search = "";
  tagsFilter: string[] = [];

  readonly headers = [
    {
      text: "Name",
      value: "name",
    },
    {
      text: "Description",
      value: "description",
    },
    {
      text: "Tags",
      value: "tags",
      filterable: false,
      sortable: false,
    },
    {
      text: "User",
      value: "member.user.name",
      filterable: false,
    },
    {
      text: "Date",
      value: "created_at",
      filterable: false,
      width: "240",
    },
    {
      text: "Actions",
      value: "action",
      sortable: false,
      filterable: false,
      width: "70",
    },
  ];

  get activeGroupId() {
    return this.groupContext.getters.activeGroupId;
  }

  get isGroupAdmin() {
    return this.groupContext.getters.isGroupAdmin;
  }

  get tags(): any[] {
    const list = this.projectsContext.getters.tags;
    return list.map((item) => {
      return {
        text: item,
      };
    });
  }

  get userProfile() {
    return this.mainContext.getters.userProfile;
  }

  get items() {
    let items = this.projectsContext.getters.projects;
    if (this.tagsFilter.length > 0) {
      items = items.filter((project) => {
        if (project.tags && this.tagsFilter.some((r) => project.tags.includes(r))) {
          return project;
        }
      });
    }
    return items;
  }

  filter(value, search, item) {
    if (!search) {
      return true;
    }
    const normalizedSearchTerm = search.toLowerCase().trim();
    if (item.name.toLowerCase().indexOf(normalizedSearchTerm) !== -1) {
      return true;
    }
    if (item.description.toLowerCase().indexOf(normalizedSearchTerm) !== -1) {
      return true;
    }
    return item.member.user.name.toLowerCase().indexOf(normalizedSearchTerm) !== -1;
  }

  removeTagFilter(item) {
    this.tagsFilter.splice(this.tagsFilter.indexOf(item), 1);
    this.tagsFilter = [...this.tagsFilter];
  }

  async deleteProject(id: number) {
    if (self.confirm("Are you sure you want to delete the project?")) {
      await this.projectsContext.actions.deleteProject(id);
    }
  }

  async mounted() {
    await Promise.all([
      this.projectsContext.actions.getGroupProjects(this.activeGroupId!),
      this.projectsContext.actions.getProjectsTags(this.activeGroupId!),
    ]);
  }
}
</script>

<style scoped>
.toolbar {
  margin-bottom: 10px;
}
</style>
